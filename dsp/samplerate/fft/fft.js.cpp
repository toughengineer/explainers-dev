#include "tools/kiss_fftr.h"

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <emscripten.h>

#include <optional>
#include <cmath>
#include <utility>
#include <limits>

#if 0
#include <iostream>
#define LOG(...) do { std::cout << __VA_ARGS__; std::cout << '\n'; } while(false)
#else
#define LOG(...) (void)0
#endif

using Float = double;

constexpr Float operator"" _float(long double x) {
  return static_cast<Float>(x);
}

enum class FftDirection {
  Direct,
  Inverse
};

struct KissFftBase {
  KissFftBase(size_t size, FftDirection direction) {
    state = kiss_fftr_alloc(size, direction == FftDirection::Direct ? 0 : 1, nullptr, nullptr);
  }
  ~KissFftBase() {
    kiss_fft_free(state);
  }
  KissFftBase(const KissFftBase &) = delete;
  KissFftBase &operator=(const KissFftBase &) = delete;

  kiss_fftr_state *state;
};

template<typename T>
struct BasicSpan {
  BasicSpan(T *data, size_t size) : data{ data }, size{ size }
  {}

  T *data;
  size_t size;
};
using Span = BasicSpan<Float>;
using ConstSpan = BasicSpan<const Float>;

Span toSpan(std::vector<Float> &data) {
  return { data.data(), data.size() };
}
ConstSpan toSpan(const std::vector<Float> &data) {
  return { data.data(), data.size() };
}

struct KissFftReal : KissFftBase {
  KissFftReal(size_t size) : KissFftBase{ size, FftDirection::Direct } {
  }

  void transform(ConstSpan input, Span output) const {
    kiss_fftr(state, input.data, reinterpret_cast<kiss_fft_cpx *>(output.data));
  }
};

struct KissFftRealBinding : KissFftReal {
  KissFftRealBinding(size_t size) : KissFftReal{ size } {
    input.resize(size);
    output.resize(size + 2);
  }

  auto transform() const {
    KissFftReal::transform(toSpan(input), toSpan(output));

    return emscripten::val{ emscripten::typed_memory_view(output.size(), output.data()) };
  }

  auto getInputTimeDataBuffer() {
    return emscripten::val{ emscripten::typed_memory_view(input.size(), input.data()) };
  }

private:
  std::vector<Float> input;
  mutable std::vector<Float> output;
};

struct KissFftRealInverse : KissFftBase {
  KissFftRealInverse(size_t size) : KissFftBase{ size, FftDirection::Inverse }
  {}

  auto transform(ConstSpan input, Span output) const {
    kiss_fftri(state, reinterpret_cast<const kiss_fft_cpx *>(input.data), output.data);
  }
};

struct KissFftRealInverseBinding : KissFftRealInverse {
  KissFftRealInverseBinding(size_t size) : KissFftRealInverse{ size } {
    input.resize(size + 2);
    output.resize(size);
  }

  auto transform() const {
    KissFftRealInverse::transform(toSpan(input), toSpan(output));

    return emscripten::val{ emscripten::typed_memory_view(output.size(), output.data()) };
  }

  auto getInputFrequencyDataBuffer() {
    return emscripten::val{ emscripten::typed_memory_view(input.size(), input.data()) };
  }

private:
  std::vector<Float> input;
  mutable std::vector<Float> output;
};


EM_JS(void, throwJSError, (const char *msg), {
  throw new Error(UTF8ToString(msg));
      });


constexpr auto pi = 3.141592653589793_float;

struct Sin {
  Sin(Float step, Float phase) :
    sinPhase{ std::sin(phase) },
    cosPhase{ std::cos(phase) },
    step{ step }
  {}
  Float next() {
    switch (state) {
    case 0:
      state = 1;
      return sinPhase;

    case 1:
      {
        state = 2;
        sinNMinus1 = std::sin(step);
        step = std::cos(step);
        cosNMinus1 = step;
        const auto result = sinNMinus1 * cosPhase + cosNMinus1 * sinPhase;
        step *= 2._float;
        return result;
      }
    }
    const auto sin = step * sinNMinus1 - sinNMinus2;
    const auto cos = step * cosNMinus1 - cosNMinus2;
    const auto result = sin * cosPhase + cos * sinPhase;
    sinNMinus2 = std::exchange(sinNMinus1, sin);
    cosNMinus2 = std::exchange(cosNMinus1, cos);
    return result;
  }

private:
  const Float sinPhase, cosPhase;
  Float
    step,
    sinNMinus2 = 0._float,
    cosNMinus2 = 1._float,
    sinNMinus1,
    cosNMinus1;
  unsigned state = 0;
};

struct Cos : Sin {
  Cos(Float step, Float phase) : Sin{ step, phase + pi / 2._float }
  {}
};

void initWindow(std::vector<Float> &window) {
  // Blackman window
  constexpr auto a = 0.16_float;
  constexpr auto a0 = 0.5_float * (1._float - a);
  constexpr auto a1 = 0.5_float;
  constexpr auto a2 = 0.5_float * a;
  const auto piDivN = pi / (window.size() + 2);
  Cos cos{ piDivN, piDivN };
  Cos cos2{ piDivN * 2._float, piDivN * 2._float };
  for (size_t n = 0; n != window.size(); ++n) {
    window[n] = a0 - a1 * cos.next() + a2 * cos2.next();
  }
}

Float approximate_log2(Float x) {
  int e;
  const auto f = std::frexp(x, &e);
  return (((-1.229145485206304_float * f + 4.824805744596478_float) * f - 7.989878313051768_float) * f + 7.840733917874801_float) * f - 3.447836997123748_float + e;
}

constexpr Float log10Of2 = 0.30102999566398114_float;

Float approximate_log10(Float x) {
  if (x == 0)
    return -std::numeric_limits<Float>::infinity();
  return approximate_log2(x) * log10Of2;
}

Float sqr(Float a) {
  return a * a;
}

Float clamp(Float value, Float low, Float high) {
  return
    value < low ? low :
    value > high ? high :
    value;
}

Float chebyshevPolynomialsDistortion(Float value, const std::vector<Float> &harmonicsCoefficients) {
  Float
    x = 1._float,
    y = value,
    sum = 0;
  if (harmonicsCoefficients.size() > 0)
    sum += harmonicsCoefficients[0] * x;
  if (harmonicsCoefficients.size() > 1)
    sum += harmonicsCoefficients[1] * y;
  if (harmonicsCoefficients.size() > 2)
    for (size_t i = 2; i != harmonicsCoefficients.size(); ++i) {
      const Float z = 2._float * value * y - x;
      sum += harmonicsCoefficients[i] * z;
      x = std::exchange(y, z);
    }
  return sum;
}

struct LowPassBiQuadFilter {
  explicit LowPassBiQuadFilter(Float normalizedFrequency) {
    const Float w0 = 2 * pi * normalizedFrequency;
    const Float cosw0 = std::cos(w0);
    const Float sinw0 = std::sin(w0);
    const Float alpha = sinw0 / 1.4142135623730951_float;  // sinw0/(2*Q)
    const Float invA0 = 1 / (1 + alpha);  // norm to a0 == 1
    a1 = -2 * cosw0 * invA0;
    a2 = (1 - alpha) * invA0;
    b0 = (1 - cosw0) / 2 * invA0;
    b1 = (1 - cosw0) * invA0;
    b2 = (1 - cosw0) / 2 * invA0;
  }
  Float operator()(Float x) {
    const auto y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
    x2 = std::exchange(x1, x);
    y2 = std::exchange(y1, y);
    return y;
  }

private:
  Float
    a1, a2, b0, b1, b2,
    x2 = 0._float, x1 = 0._float, y2 = 0._float, y1 = 0._float;
};

struct Oversampling {
  Oversampling()
  {}

  size_t getSampleRate() const {
    return sampleRate;
  }
  void setSampleRate(size_t rate) {
    if (sampleRate != rate) {
      sampleRate = rate;
    }
  }

  size_t getOversample() const {
    return oversample;
  }
  void setOversample(size_t multiplier) {
    if (multiplier == 0)
      throwJSError("oversample value must be >= 1");

    if (oversample == multiplier)
      return;

    oversample = multiplier;
    if (oversample == 1)
    {
      oversampledFftInput = {};
      oversampledFftOutput = {};
      oversampledFrequencyData = {};
      oversampledWindow = {};
    }
  }

  Float getCutoffFrequency() const {
    return cutoffFrequency;
  }
  void setCutoffFrequency(Float frequency) {
    cutoffFrequency = frequency;
  }

  size_t getFftSize() const {
    return fftSize;
  }
  void setFftSize(size_t size) {
    if (size == 0 || (size & 1) == 1)
      throwJSError("FFT size must be greater than 0 and even");

    if (fftSize != size) {
      fftSize = size;
      if (fft)
        fft.reset();
    }
  }

  Float getStartFrequency() const {
    return startFrequency;
  }
  void setStartFrequency(Float f) {
    if (f <= 0._float)
      throwJSError("Frequency must be greater than 0");

    startFrequency = f;
  }

  Float getEndFrequency() const {
    return startFrequency;
  }
  void setEndFrequency(Float f) {
    if (f <= 0._float)
      throwJSError("Frequency must be greater than 0");

    endFrequency = f;
  }

  auto getHarmonicsCoefficients() const {
    return emscripten::val::array(harmonicsCoefficients);
  }
  void setHarmonicsCoefficients(emscripten::val coefficients) {
    harmonicsCoefficients = emscripten::convertJSArrayToNumberVector<Float>(coefficients);
  }

  auto getFrequencyData() const {
    return emscripten::val{ emscripten::typed_memory_view(frequencyData.size(), frequencyData.data()) };
  }
  auto getOversampledFrequencyData() const {
    return emscripten::val{ oversample == 1 ? emscripten::typed_memory_view(frequencyData.size(), frequencyData.data()) :
      emscripten::typed_memory_view(oversampledFrequencyData.size(), oversampledFrequencyData.data()) };
  }

  Float getMinDecibels() const {
    return dBMin;
  }
  void setMinDecibels(Float dB) {
    dBMin = dB;
  }

  Float getMaxDecibels() const {
    return dBMax;
  }
  void setMaxDecibels(Float dB) {
    dBMax = dB;
  }

  auto getAudioBuffer() const {
    return emscripten::val{ emscripten::typed_memory_view(audioBuffer.size(), audioBuffer.data()) };
  }

  size_t getAudioBufferSampleCount() const {
    return audioBufferSampleCount;
  }

  void restartAndRenderStep() {
    chunkSize = sampleRate * lengthInSeconds / numberOfChunks;

    const auto oversampledChunkSize = oversample * chunkSize;

    const auto oversampledTimeStep = 1._float / (oversample * sampleRate);

    const auto safeEndFrequency = std::min(endFrequency, sampleRate / 2._float);

    // x(t) = sin(2*pi*(f0*t + (f1 - f0)/(2*T)*t^2))
    u1 = 2._float * pi * startFrequency * oversampledTimeStep;
    u2 = pi * (safeEndFrequency - startFrequency) / (oversampledChunkSize * numberOfChunks) * oversampledTimeStep;

    LOG("restarting with\n"
        "Fs: " << sampleRate << " Hz\n"
        "sweep start: " << startFrequency << " Hz\n"
        "      end: " << safeEndFrequency << " Hz\n"
        "fft size: " << fftSize << "\n"
        "oversample: " << oversample << "\n"
        "cutoff frequency: " << cutoffFrequency << " Hz\n"
        "chunk size: " << chunkSize << "\n"
        "duration: " << (oversampledChunkSize * numberOfChunks * oversampledTimeStep) << "\n"
    );
    LOG("harmonics coefficients (" << harmonicsCoefficients.size() << "):");
    for (auto c : harmonicsCoefficients)
      LOG("  " << c);

    if (2 * oversampledChunkSize > fftSize)
      throwJSError("FFT size is too small");

    if (!fft)
    {
      fft.emplace(fftSize);
      fftNorm = 1._float / fftSize;
    }

    auto initFftInput = [](auto &input, size_t chunkSize) {
      for (size_t i = 0; i != chunkSize; ++i)
        input[i] = 0._float;
      for (size_t i = 2 * chunkSize; i != input.size(); ++i)
        input[i] = 0._float;
    };

    fftInput.resize(fftSize);
    initFftInput(fftInput, chunkSize);

    fftOutput.resize(fftSize + 2);
    frequencyData.resize(fftSize / 2);
    window.resize(chunkSize - 1);
    initWindow(window);
    if (oversample != 1)
    {
      oversampledFftInput.resize(fftSize);
      initFftInput(oversampledFftInput, oversampledChunkSize);

      oversampledFftOutput.resize(fftSize + 2);
      oversampledFrequencyData.resize(fftSize / 2);
      oversampledWindow.resize(oversampledChunkSize - 1);
      initWindow(oversampledWindow);
    }

    filter = LowPassBiQuadFilter(cutoffFrequency * oversampledTimeStep);
    filter2 = LowPassBiQuadFilter(cutoffFrequency * oversampledTimeStep);

    audioBuffer.clear();
    audioBufferSampleCount = chunkSize * numberOfChunks;
    audioBuffer.reserve(audioBufferSampleCount);

    chunkBuffer.resize(oversampledChunkSize, 0._float);

    dBDelta = 255._float / (dBMax - dBMin);

    if (!harmonicsCoefficients.empty()) {
      harmonicsCoefficientsNorm = 0._float;
      for (auto coefficient : harmonicsCoefficients)
        harmonicsCoefficientsNorm += coefficient;
      harmonicsCoefficientsNorm = 1._float / harmonicsCoefficientsNorm;
    }

    renderStep();
  }

  void renderNextStep() {
    if (oversample == 1) {
      for (size_t i = 0; i != window.size(); ++i)
        fftInput[i] = chunkBuffer[i] * window[i];
      fftInput[chunkSize - 1] = chunkBuffer[chunkSize - 1];
    }
    else {
      {
        for (size_t i = 0; i != window.size(); ++i)
          fftInput[i] = chunkBuffer[i * oversample] * window[i];
        fftInput[chunkSize - 1] = chunkBuffer[(chunkSize - 1) * oversample];
      }
      for (size_t i = 0; i != oversampledWindow.size(); ++i)
        oversampledFftInput[i] = chunkBuffer[i] * oversampledWindow[i];
      oversampledFftInput[oversampledWindow.size()] = chunkBuffer[oversampledWindow.size()];
    }

    renderStep();
  }

private:
  void renderStep() {
    const auto oversampledChunkSize = oversample * chunkSize;

    Float currentSampleIndex = audioBuffer.size() * oversample;
    for (size_t i = 0; i != oversampledChunkSize; ++i) {
      auto sample = std::sin((u1 + u2 * currentSampleIndex) * currentSampleIndex);
      if (!harmonicsCoefficients.empty())
        sample = chebyshevPolynomialsDistortion(sample, harmonicsCoefficients) * harmonicsCoefficientsNorm;
      chunkBuffer[i] = clamp(filter2(filter(sample)), -1._float, 1._float);
      currentSampleIndex += 1._float;
    }

    if (oversample == 1) {
      audioBuffer.insert(audioBuffer.end(), chunkBuffer.begin(), chunkBuffer.end());
      fftInput[chunkSize] = chunkBuffer[0];
      for (size_t i = 1; i != chunkSize; ++i) {
        fftInput[chunkSize + i] = chunkBuffer[i] * window[window.size() - i];
      }
      estimatePowerSpectrumDensity(fftInput, fftOutput, frequencyData);
    }
    else {
      {
        {
          const auto sample = chunkBuffer[0];
          audioBuffer.push_back(sample);
          fftInput[chunkSize] = sample;
        }
        for (size_t i = 1; i != chunkSize; ++i) {
          const auto sample = chunkBuffer[i * oversample];
          audioBuffer.push_back(sample);
          fftInput[chunkSize + i] = sample * window[window.size() - i];
        }
        estimatePowerSpectrumDensity(fftInput, fftOutput, frequencyData);
      }
      oversampledFftInput[oversampledChunkSize] = chunkBuffer[0];
      for (size_t i = 1; i != oversampledChunkSize; ++i)
        oversampledFftInput[oversampledChunkSize + i] = chunkBuffer[i] * oversampledWindow[oversampledWindow.size() - i];
      estimatePowerSpectrumDensity(oversampledFftInput, oversampledFftOutput, oversampledFrequencyData);
    }
  }

  void estimatePowerSpectrumDensity(const std::vector<Float> &input,
                                    std::vector<Float> &output,
                                    std::vector<uint8_t> &powerSpectrumDensity) {
    fft->transform(toSpan(input), toSpan(output));
    auto dBToUint8 = [this](Float dB) {
      const auto dBMapped = dBDelta * (dB - dBMin);
      return
        dBMapped < 0._float ? 0 :
        dBMapped > 255._float ? 255 :
        static_cast<uint8_t>(std::round(dBMapped));
    };
    powerSpectrumDensity[0] = dBToUint8(20._float * approximate_log10(fftNorm * std::abs(output[0])));
    const auto sqrNorm = fftNorm * fftNorm;
    for (size_t i = 1; i != powerSpectrumDensity.size(); ++i) {
      const auto sqrSum = (sqr(output[i * 2]) + sqr(output[i * 2 + 1]));
      const auto dB = 20._float * (0.5_float * approximate_log10(sqrNorm * sqrSum) + log10Of2);  // 20*log10(2*sqrSum^0.5/N)
      powerSpectrumDensity[i] = dBToUint8(dB);
    }
  }

  static constexpr auto lengthInSeconds = 2;
  static constexpr auto numberOfChunks = 1000;
  static constexpr auto defaultSampleRate = 48000;

  std::vector<Float>
    fftInput,
    fftOutput,
    oversampledFftInput,
    oversampledFftOutput,
    chunkBuffer,
    window,
    oversampledWindow,
    harmonicsCoefficients;
  std::vector<uint8_t>
    frequencyData,
    oversampledFrequencyData;
  std::vector<float> audioBuffer;
  std::optional<KissFftReal> fft;
  LowPassBiQuadFilter filter{ defaultSampleRate / 2._float };
  LowPassBiQuadFilter filter2{ defaultSampleRate / 2._float };
  Float
    cutoffFrequency = defaultSampleRate / 2._float,
    startFrequency = 55._float,
    endFrequency = 14080._float,
    dBMin = -100._float,
    dBMax = -30._float,
    dBDelta,
    harmonicsCoefficientsNorm,
    u1, u2,
    fftNorm;
  size_t
    fftSize = 2048,
    sampleRate = defaultSampleRate,
    oversample = 1,
    chunkSize,
    audioBufferSampleCount;
};

EMSCRIPTEN_BINDINGS(KissFft) {
  emscripten::class_<KissFftRealBinding>("KissFftReal")
    .constructor<size_t>()
    .function("getInputTimeDataBuffer", &KissFftRealBinding::getInputTimeDataBuffer)
    .function("transform", &KissFftRealBinding::transform)
    ;

  emscripten::class_<KissFftRealInverseBinding>("KissFftRealInverse")
    .constructor<size_t>()
    .function("getInputFrequencyDataBuffer", &KissFftRealInverseBinding::getInputFrequencyDataBuffer)
    .function("transform", &KissFftRealInverseBinding::transform)
    ;

  emscripten::class_<Oversampling>("Oversampling")
    .constructor()
    .property("sampleRate", &Oversampling::getSampleRate, &Oversampling::setSampleRate)
    .property("oversample", &Oversampling::getOversample, &Oversampling::setOversample)
    .property("cutoffFrequency", &Oversampling::getCutoffFrequency, &Oversampling::setCutoffFrequency)
    .property("fftSize", &Oversampling::getFftSize, &Oversampling::setFftSize)
    .property("startFrequency", &Oversampling::getStartFrequency, &Oversampling::setStartFrequency)
    .property("endFrequency", &Oversampling::getEndFrequency, &Oversampling::setEndFrequency)
    .property("harmonicsCoefficients", &Oversampling::getHarmonicsCoefficients, &Oversampling::setHarmonicsCoefficients)
    .function("restartAndRenderStep", &Oversampling::restartAndRenderStep)
    .function("renderNextStep", &Oversampling::renderNextStep)
    .property("stepFrequencyData", &Oversampling::getFrequencyData)
    .property("stepOversampledFrequencyData", &Oversampling::getOversampledFrequencyData)
    .property("minDecibels", &Oversampling::getMinDecibels, &Oversampling::setMinDecibels)
    .property("maxDecibels", &Oversampling::getMaxDecibels, &Oversampling::setMaxDecibels)
    .property("audioBufferData", &Oversampling::getAudioBuffer)
    .property("audioBufferSampleCount", &Oversampling::getAudioBufferSampleCount)
    ;
}
