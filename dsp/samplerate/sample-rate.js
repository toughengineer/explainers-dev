let debug = false;

const graphStyle = {
  line: {
    color: '#0882ff',
    color2: '#ffa631',
    width: 3
  },
  xAxis: {
    color: '#aaa'
  },
  yAxis: {
    grid: {
      color: '#ddd'
    }
  },
  mark: {
    line: {
      color: '#f20a0a',
      width: 3
    },
    fillOpacity: 0.85,
    radius: 5
  }
};

let Module = {
  preRun: [],
  postRun: [],
  print: function (text) {
    if (arguments.length > 1)
      text = Array.prototype.slice.call(arguments).join(' ');
    console.log(text);
  },
  printErr: function (text) {
    if (arguments.length > 1)
      text = Array.prototype.slice.call(arguments).join(' ');
    console.error(text);
  },
  setStatus: function (text) {
    if (!Module.setStatus.last)
      Module.setStatus.last = { time: Date.now(), text: '' };
    if (text === Module.setStatus.last.text)
      return;
    console.log(text);
  },
  totalDependencies: 0,
  monitorRunDependencies: function (left) {
    this.totalDependencies = Math.max(this.totalDependencies, left);
    Module.setStatus(left ? 'Preparing... (' + (this.totalDependencies - left) + '/' + this.totalDependencies + ')' : 'All downloads complete.');
  }
};

function* genSine(step, phase = 0.) {
  //sin(a + b) = sin(a)*cos(b) + cos(a)*sin(b)
  //sin(n*x) = 2*cos(x)*sin((n - 1)*x) - sin((n - 2)*x)
  //cos(n*x) = 2*cos(x)*cos((n - 1)*x) - cos((n - 2)*x)
  const sinPhase = Math.sin(phase);
  const cosPhase = Math.cos(phase);
  //a == 0
  let sinNMinus2 = 0.;
  let cosNMinus2 = 1.;
  yield /*sinNMinus2 * cosPhase + cosNMinus2 * */sinPhase;
  //a == 1*step
  let sinNMinus1 = Math.sin(step);
  const cosStep = Math.cos(step);
  let cosNMinus1 = cosStep;
  yield sinNMinus1 * cosPhase + cosNMinus1 * sinPhase;
  //a == n*step
  const cosStep2 = 2 * cosStep;
  while (true) {
    const sin = cosStep2 * sinNMinus1 - sinNMinus2;
    const cos = cosStep2 * cosNMinus1 - cosNMinus2;
    yield sin * cosPhase + cos * sinPhase;
    sinNMinus2 = sinNMinus1;
    sinNMinus1 = sin;
    cosNMinus2 = cosNMinus1;
    cosNMinus1 = cos;
  }
}

function getGenSine(step, phase = 0.) {
  const sine = genSine(step, phase);
  return function () {
    return sine.next().value;
  }
}

function sCurve(x) {
  x = 2 * x - 1;
  return x / Math.pow(1 + 16 * x * x * x * x, 1 / 4) + 0.5;
}

function getDrawOnNextAnimationFrameThrottled(drawCallback) {
  let requested = false;
  return () => {
    if (!requested) {
      requestAnimationFrame(() => {
        drawCallback();
        requested = false;
      });
      requested = true;
    }
  };
}

function getThrottledCallback(callback, period, delay = 0) {
  let requested = false;
  let scheduled = false;
  function request() {
    if (scheduled) {
      requested = true;
    }
    else {
      scheduled = true;
      const timeout = (() => {
        if (!requested && delay != 0) {
          requested = true;
          return delay;
        }

        callback();
        return period;
      })();
      setTimeout(() => {
        scheduled = false;
        if (requested) {
          request();
          requested = false;
        }
      }, timeout);
    }
  };
  return request;
}

(() => {
  const ctx = document.getElementById('sampling').getContext('2d');
  const parameter = document.getElementById('samplingParameter');

  function draw() {
    const param = Number.parseFloat(parameter.value) / 100;

    ctx.save();
    ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
    ctx.setTransform(1, 0, 0, -1, 0.5, Math.trunc(ctx.canvas.height / 2) + 0.5);

    ctx.strokeStyle = graphStyle.xAxis.color;
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(ctx.canvas.width - 1, 0);
    ctx.stroke();

    const step = 0.02;
    const m = ctx.canvas.height / 2;

    const sin1 = getGenSine(step / 2, -param * 5);
    const m1 = 0.9 * m * (1 - sCurve(1.25 * param));
    const sin2 = getGenSine(step, -3.5 - param * 7);
    const s3 = sCurve(2 * (param - 0.45));
    const m2 = m * (sCurve(2 * (param - 0.25)) - s3 * 0.5);
    const sin3 = getGenSine(step * 3, -3.1 - param * 22);
    const m3 = 0.5 * m * s3;
    function f() {
      return sin1() * m1 +
        sin2() * m2 +
        sin3() * m3;
    }

    let points = [];

    ctx.strokeStyle = graphStyle.line.color;
    ctx.lineWidth = graphStyle.line.width;
    ctx.lineCap = 'round';
    ctx.beginPath();
    ctx.moveTo(0, f());
    let counter = 12;
    for (let i = 1; i != ctx.canvas.width / 2 + 1; ++i) {
      const x = i * 2;
      const y = f();
      ctx.lineTo(x, y);
      if (++counter == 25) {
        points.push([x, y]);
        counter = 0;
      }
    }
    ctx.stroke();

    const radius = graphStyle.mark.radius;
    ctx.strokeStyle = graphStyle.yAxis.grid.color;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let p of points) {
      if (Math.abs(p[1]) > radius) {
        const offset = Math.sign(p[1]);
        ctx.moveTo(p[0], offset);
        ctx.lineTo(p[0], p[1] - offset * radius);
      }
    }
    ctx.stroke();

    ctx.strokeStyle = graphStyle.mark.line.color;
    ctx.lineWidth = graphStyle.mark.line.width;
    ctx.fillStyle = graphStyle.line.color2;
    for (let p of points) {
      ctx.beginPath();
      ctx.ellipse(p[0], p[1], radius, radius, 0, 0, 2 * Math.PI);
      ctx.globalAlpha = graphStyle.mark.fillOpacity;
      ctx.fill();
      ctx.globalAlpha = 1;
      ctx.stroke();
    }


    ctx.restore();
  }
  draw();

  parameter.addEventListener('input', getDrawOnNextAnimationFrameThrottled(draw));
})();


function fillCanvas(ctx, color = 'black') {
  ctx.save();

  ctx.fillStyle = color;
  ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);

  ctx.restore();
}

(() => {
  const ctx = document.getElementById('fourierSeries').getContext('2d', { alpha: false });
  const componentCount = document.getElementById('fourierSeriesComponentsCount');
  const countDisplay = document.getElementById('fourierSeriesCountDisplay');

  // step at 1 on [0..2], Period = 2
  // f(t) = A_0/2 + sum_1_inf(A_n*cos(2*pi*n*t/P - phase_n)
  // a_n = 2/P * inegral_P(f(x)*cos(2*pi*n/P * x) dx) = integral_1_2(cos(pi*n * x) dx) = (sin(2*pi*n) - sin(pi*n))/(pi*n)
  // a_0 = 1
  // a_n = 0 for n >= 1
  // b_n = 2/P * inegral_P(f(x)*sin(2*pi*n/P * x) dx) = integral_1_2(sin(pi*n * x) dx) = (cos(pi*n) - cos(2*pi*n))/(pi*n) = (cos(pi*n) - 1)/(pi*n)
  // b_0 = 0
  // b_n = -2/(pi*n) for odd n, 0 otherwise (0, -2/pi, 0, -2/(3*pi), 0, -2/(5*pi) ...)
  // A_n = sqrt(a_n^2 + b_n^2) = |b_n|
  // phase_n = atan2(b_n, a_n) = pi/2
  // f(t) = 1/2 + sum_1_inf(A_n*cos(pi*n*t - pi/2) = 1/2 + sum_1_inf(A_n*sin(pi*n*t))

  function draw() {
    const x1 = 1;

    fillCanvas(ctx, 'white');

    ctx.save();
    ctx.setTransform(1, 0, 0, -1, 0.5, Math.trunc(ctx.canvas.height * 7 / 8) + 0.5);

    ctx.strokeStyle = graphStyle.xAxis.color;
    ctx.lineWidth = 1;

    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(ctx.canvas.width - 1, 0);
    ctx.stroke();

    ctx.strokeStyle = graphStyle.line.color;
    ctx.lineWidth = 7;

    const halfWidth = ctx.canvas.width / 2;
    const yScale = ctx.canvas.height * 6 / 8;

    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(halfWidth, 0);
    ctx.lineTo(halfWidth, yScale);
    ctx.lineTo(ctx.canvas.width - 1, yScale);
    ctx.stroke();

    ctx.strokeStyle = '#333';
    ctx.lineWidth = 1;
    ctx.lineCap = 'round';
    ctx.lineJoin = 'round';

    const values = new Array(Math.round(ctx.canvas.width) + 1).fill(0.5);
    const a = -2 / Math.PI;
    function calcNDrawNext(n, skip) {
      const A = a / n;
      const k = Math.PI * n;
      const sin = getGenSine(k * x1 / ctx.canvas.width, k * 0.5);

      ctx.beginPath();
      {
        const y = sin() * A;
        values[0] += y;
        ctx.moveTo(0, (0.5 + y) * yScale);
      }
      skip += 1;
      for (let i = 1; i != values.length - 1; ++i) {
        const y = sin() * A;
        values[i] += y;
        if (skip < 2 || i % skip == 0)
          ctx.lineTo(i, (0.5 + y) * yScale);
      }
      {
        const y = sin() * A;
        values[values.length - 1] += y;
        ctx.lineTo(values.length - 1, (0.5 + y) * yScale);
      }
      ctx.stroke();
    }
    function drawComponent() {
      ctx.beginPath();
      ctx.moveTo(0, values[0] * yScale);
      for (let i = 1; i != values.length; ++i) {
        ctx.lineTo(i, values[i] * yScale);
      }
      ctx.stroke();
    }

    const N = Number.parseInt(componentCount.value);
    for (let m = 1; m != N; ++m) {
      ctx.globalAlpha = 0.9 / (N - m + 1) + 0.1;
      calcNDrawNext(2 * m - 1, Math.round(m / 100));
    }

    ctx.strokeStyle = graphStyle.line.color2;
    calcNDrawNext(2 * N - 1, 0);
    ctx.globalAlpha = 1;
    ctx.lineWidth = graphStyle.line.width;
    drawComponent();

    ctx.restore();

    countDisplay.textContent = N.toString();
  }
  draw();

  componentCount.addEventListener('input', getDrawOnNextAnimationFrameThrottled(draw));
})();


function lerp(a, b, t) {
  return a * (1 - t) + b * t;
}

const lanczosKernel = (() => {
  const a = 5;
  const coefficients = new Array(146);
  const indexNorm = a / (coefficients.length - 1);
  coefficients[0] = 1;
  coefficients[coefficients.length - 1] = 0;
  for (let i = 1; i != coefficients.length - 1; ++i) {
    const z = i * indexNorm * Math.PI;
    coefficients[i] = a * Math.sin(z) * Math.sin(z / a) / (z * z);
  }

  const valueNorm = (coefficients.length - 1) / a;
  return (x) => {
    x = Math.abs(x);
    if (x >= a)
      return 0;
    const y = x * valueNorm;
    const i = Math.floor(y);
    const f = y - i;
    return lerp(coefficients[i], coefficients[i + 1], f);
  };
})();


function rand(x) {
  function lcgNext(x) {
    const a = 1664525;
    const c = 1013904223;
    // m = 2^32
    return ((a * x) & 0xffffffff + c) & 0xffffffff;
  }
  return lcgNext(lcgNext(x)) / 0x7fffffff * 4 - 1;
}

function noise(x) {
  const i = Math.floor(x);
  const f = x - i;
  if (debug)
    return lerp(rand(i), rand(i + 1), f);
  let s = 0;
  for (let k = -5; k != 6; ++k) {
    s += lanczosKernel(k - f) * rand(i + k);
  }
  return s;
}

function blackmanWindow(x) {
  // alpha = 0.16;
  const a0 = 0.42;
  const a1 = 0.5;
  const a2 = 0.08;
  const y = 2 * Math.PI * x;
  return a0 - a1 * Math.cos(y) + a2 * Math.cos(2 * y);
}

(() => {
  const totalNumberOfSamples = 64;
  const numberOfSamples = 20;
  const numberOfSamplesBefore = Math.floor((totalNumberOfSamples - numberOfSamples) / 2);
  const numberOfSamplesAfter = totalNumberOfSamples - numberOfSamples - numberOfSamplesBefore;

  const shapePicker = document.getElementById('waveShape');
  const ctx = document.getElementById('reconstruction').getContext('2d', { alpha: false });
  const frequencySlider = document.getElementById('waveFrequency');
  const randomizeNoiseButton = document.getElementById('randomizeNoiseButton');

  const yScale = 0.75 * ctx.canvas.height / 2;
  const sampleDistance = Math.round(ctx.canvas.width / numberOfSamples);

  const stepsPerPixel = 2;

  const shapes = new Map();
  shapes.set('sine', frequency => {
    const k = 2 * Math.PI * frequency / ctx.canvas.width;
    const drawingStep = k / stepsPerPixel;
    return {
      func: x => Math.sin(k * x),
      drawing: getGenSine(drawingStep)
    }
  });

  function getFuncSteps(f, frequency) {
    let x = 0;
    const step = frequency / ctx.canvas.width;
    const drawingStep = step / stepsPerPixel;
    return {
      func: x => f(step * x),
      drawing: () => {
        const y = f(x);
        x += drawingStep;
        return y;
      }
    };
  }

  const mod1 = x => x == 0 ? -1 : x % 1;
  const sign = x => x == 0 ? -1 : Math.sign(x);
  shapes.set('square', frequency => getFuncSteps(x => sign(sign(x) * 0.5 - mod1(x)), frequency));
  shapes.set('triangle', frequency => getFuncSteps(x => 1 - 4 * Math.abs(mod1(x + 0.25) - sign(x + 0.25) * 0.5), frequency));
  shapes.set('saw', frequency => getFuncSteps(x => 2 * mod1(x) - sign(x), frequency));
  shapes.set('saw2', frequency => getFuncSteps(x => sign(x) - 2 * mod1(x), frequency));
  let noiseSeed = 5;
  shapes.set('noise', frequency => getFuncSteps(x => noise(2 * x + noiseSeed), frequency));

  let drawReconstructed = (ctx) => {
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.fillStyle = 'gray';
    ctx.textBaseline = 'top';
    ctx.textAlign = 'right';
    ctx.font = '16px sans-serif';
    ctx.fillText('loading module...', ctx.canvas.width - 5, 5);
  };

  function draw() {

    fillCanvas(ctx, 'white');

    ctx.save();
    ctx.setTransform(1, 0, 0, -1, 0.5, Math.trunc(ctx.canvas.height / 2) + 0.5);

    ctx.strokeStyle = graphStyle.xAxis.color;
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(ctx.canvas.width - 1, 0);
    ctx.stroke();

    try {

      const f = (shapes.get(shapePicker.value))(Number.parseFloat(frequencySlider.value));

      const samples = [];
      const xs = [];

      const windowSize = 0.66;
      for (let i = -numberOfSamplesBefore; i != 0; ++i) {
        const x = i * sampleDistance;
        xs.push(x);
        let sample = f.func(x);
        const z = 1 + i / numberOfSamplesBefore;
        if (z < windowSize)
          sample *= blackmanWindow(0.5 * z / windowSize);
        samples.push(sample);
      }

      ctx.strokeStyle = graphStyle.line.color;
      ctx.lineCap = 'round';
      ctx.lineJoin = 'round';

      ctx.beginPath();
      const y = f.drawing();
      xs.push(0);
      samples.push(y);
      ctx.moveTo(0, y * yScale);
      ctx.lineTo(0.5, f.drawing() * yScale);
      for (let x = 1; x != ctx.canvas.width; ++x) {
        const y = f.drawing();
        if (x % sampleDistance == 0) {
          xs.push(x);
          samples.push(y);
        }
        ctx.lineTo(x, y * yScale);
        ctx.lineTo(x + 0.5, f.drawing() * yScale);
      }
      ctx.lineWidth = graphStyle.line.width * 3;
      ctx.globalAlpha = 0.2;
      ctx.stroke();
      ctx.lineWidth = 1;
      ctx.globalAlpha = 1;
      ctx.stroke();

      const radius = graphStyle.mark.radius;
      ctx.strokeStyle = graphStyle.yAxis.grid.color;
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (let i = numberOfSamplesBefore + 1; i != samples.length; ++i) {
        if (xs[i] >= ctx.canvas.width)
          break;
        const y = samples[i] * yScale;
        if (Math.abs(y) > radius) {
          const offset = Math.sign(y);
          ctx.moveTo(xs[i], offset);
          ctx.lineTo(xs[i], y - offset * radius);
        }
      }
      ctx.stroke();

      ctx.strokeStyle = graphStyle.mark.line.color;
      ctx.lineWidth = graphStyle.mark.line.width;
      ctx.fillStyle = graphStyle.line.color2;
      for (let i = numberOfSamplesBefore + 1; i != samples.length; ++i) {
        if (xs[i] >= ctx.canvas.width)
          break;
        ctx.beginPath();
        ctx.ellipse(xs[i], samples[i] * yScale, radius, radius, 0, 0, 2 * Math.PI);
        ctx.globalAlpha = graphStyle.mark.fillOpacity;
        ctx.fill();
        ctx.globalAlpha = 1;
        ctx.stroke();
      }

      const endX = xs[xs.length - 1];
      for (let i = 0; i != numberOfSamplesAfter; ++i) {
        const x = endX + (i + 1) * sampleDistance;
        xs.push(x);
        let sample = f.func(x);
        const z = 1 - i / numberOfSamplesAfter;
        if (z < windowSize)
          sample *= blackmanWindow(0.5 * z / windowSize);
        samples.push(sample);
      }

      drawReconstructed(ctx, samples);

    }
    catch (e) {
      console.error(`Error occurred while drawing shape '${shapePicker.selectedOptions[0].textContent}' (${shapePicker.value}): ${e}`);
    }
    ctx.restore();
  }
  draw();

  Module.postRun.push(() => {
    const numberOfUpsampledSamples = 1024;
    const fft = new Module.KissFftReal(totalNumberOfSamples);
    const ifft = new Module.KissFftRealInverse(numberOfUpsampledSamples);
    drawReconstructed = (ctx, samples) => {
      const input = fft.getInputTimeDataBuffer();
      input.set(samples);
      const frequencies = fft.transform();
      const invInput = ifft.getInputFrequencyDataBuffer();
      for (let i = 0; i != totalNumberOfSamples; ++i)
        invInput[i] = frequencies[i];
      invInput[totalNumberOfSamples] = 0.5 * frequencies[totalNumberOfSamples];
      invInput[totalNumberOfSamples + 1] = 0.5 * frequencies[totalNumberOfSamples + 1];

      const upsampled = ifft.transform();

      const scale = yScale / totalNumberOfSamples;

      ctx.strokeStyle = graphStyle.line.color2;
      ctx.lineWidth = graphStyle.line.width;

      ctx.beginPath();
      const xScale =
        debug ? numberOfUpsampledSamples / ctx.canvas.width :
          ctx.canvas.width * totalNumberOfSamples / (numberOfSamples * numberOfUpsampledSamples);
      const first = debug ? 0 : numberOfUpsampledSamples * numberOfSamplesBefore / totalNumberOfSamples;
      let i = first;
      ctx.moveTo(0, upsampled[i] * scale);
      for (++i; i != upsampled.length; ++i) {
        const x = (i - first) * xScale;
        ctx.lineTo(x, upsampled[i] * scale);
        if (x > ctx.canvas.width)
          break;
      }
      ctx.stroke();
    };
    draw();
  });



  frequencySlider.addEventListener('input', getDrawOnNextAnimationFrameThrottled(draw));

  shapePicker.addEventListener('change', () => {
    draw();
    if (shapePicker.value == 'noise')
      randomizeNoiseButton.classList.remove('hidden');
    else
      randomizeNoiseButton.classList.add('hidden');
  });

  randomizeNoiseButton.addEventListener('click', () => {
    noiseSeed = (Math.random() * 0xffffffff) & 0xffffffff;
    draw();
  });
})();


const getDefaultSampleRate = (() => {
  try {
    const audioCtx = new AudioContext();
    const sampleRate = audioCtx.sampleRate;
    audioCtx.close();
    return function () { return sampleRate; }
  }
  catch (e) {
    return function () { return null; }
  }
})();

function setDefaulSampleRateOption(e) {
  const defaultSampleRate = getDefaultSampleRate();

  let option = null;
  if (defaultSampleRate !== null)
    option = e.querySelector(`option[value="${defaultSampleRate}"]`);

  if (option !== null) {
    option.textContent += ' (default)';
  }
  else {
    option = document.createElement('OPTION');
    if (defaultSampleRate !== null) {
      option.value = defaultSampleRate.toString();
      option.textContent = `default (${defaultSampleRate} Hz)`;
    }
    else {
      option.value = '';
      option.textContent = `default`;
    }
    e.insertBefore(option, e.childNodes[0]);
  }
  option.selected = true;
}

function clamp(value, min, max) {
  return value < min ? min : value > max ? max : value;
}

function setupColorSchemeChooser(rootSelector, ctx, callbackOnChange) {
  const colorOption = document.querySelector(rootSelector + ' input.color');
  const whiteOnBlackOption = document.querySelector(rootSelector + ' input.white');
  const blackOnWhiteOption = document.querySelector(rootSelector + ' input.black');

  function getIntensityColor(intensity) {
    if (intensity == 255)
      return [255, 255, 255];
    const r = clamp((intensity - 64) * 4, 0, 255);
    const g = Math.max((intensity - 192) * 4, 0);
    const b = clamp(383 - Math.abs(intensity - 96) * 4, 0, 255);
    return [r, g, b];
  }

  function getIntensityWhiteOnBlack(intensity) {
    return [intensity, intensity, intensity];
  }

  function getIntensityBlackOnWhite(intensity) {
    const value = 255 - intensity;
    return [value, value, value];
  }

  const colorSpectrumGradient = (() => {
    const g = ctx.createLinearGradient(0, 0, 0, ctx.canvas.height);
    g.addColorStop(0, 'rgb(63, 63, 63)');
    g.addColorStop(0.25, 'rgb(0, 0, 255)');
    g.addColorStop(0.5, 'rgb(255, 0, 255)');
    g.addColorStop(0.75, 'rgb(255, 0, 0)');
    g.addColorStop(0.99, 'rgb(255, 255, 0)');
    g.addColorStop(1, 'rgb(255, 255, 255)');
    return g;
  })();

  const whiteOnBlackSpectrumGradient = (() => {
    const g = ctx.createLinearGradient(0, 0, 0, ctx.canvas.height);
    g.addColorStop(0, 'rgb(63, 63, 63)');
    g.addColorStop(0.25, 'rgb(63, 63, 63)');
    g.addColorStop(1, 'rgb(255, 255, 255)');
    return g;
  })();

  const blackOnWhiteSpectrumGradient = (() => {
    const g = ctx.createLinearGradient(0, 0, 0, ctx.canvas.height);
    g.addColorStop(0, 'rgb(193, 193, 193)');
    g.addColorStop(0.25, 'rgb(193, 193, 193)');
    g.addColorStop(1, 'rgb(0, 0, 0)');
    return g;
  })();

  function getScheme() {
    if (whiteOnBlackOption.checked)
      return {
        intensity: getIntensityWhiteOnBlack,
        gradient: whiteOnBlackSpectrumGradient,
        baseColor: 'black'
      };
    else if (blackOnWhiteOption.checked)
      return {
        intensity: getIntensityBlackOnWhite,
        gradient: blackOnWhiteSpectrumGradient,
        baseColor: 'white'
      };
    // else
    return {
      intensity: getIntensityColor,
      gradient: colorSpectrumGradient,
      baseColor: 'black'
    };
  }

  const handler = () => {
    callbackOnChange(getScheme());
  };

  colorOption.addEventListener('change', handler);
  whiteOnBlackOption.addEventListener('change', handler);
  blackOnWhiteOption.addEventListener('change', handler);

  return getScheme();
}

function formatHz(value) {
  if (value % 100 == 0)
    return `${value / 1000} kHz`;
  return `${value} Hz`;
}

class LinearSpectrogram {
  constructor(ctx) {
    this.slice = ctx.createImageData(1, ctx.canvas.height);
  }
  drawSlice(ctx, x, frequencyData, intensityToColor) {
    ctx.save();
    const height = ctx.canvas.height;
    const data = this.slice.data;
    const pixelsPerSamples = height / frequencyData.length;
    let y = height - 1;
    let currentSum = 0;
    let currentIntensity = 0;
    for (let frequencySample of frequencyData) {
      currentSum += pixelsPerSamples;
      currentIntensity = Math.max(currentIntensity, frequencySample);
      if (currentSum >= 1) {
        const [r, g, b] = intensityToColor(currentIntensity);
        data[y * 4] = r;
        data[y * 4 + 1] = g;
        data[y * 4 + 2] = b;
        data[y * 4 + 3] = 255;
        --y;
        currentSum -= 1;
        currentIntensity = 0;
      }
    }
    if (currentSum >= 1) {
      const [r, g, b] = intensityToColor(currentIntensity);
      data[y * 4] = r;
      data[y * 4 + 1] = g;
      data[y * 4 + 2] = b;
      data[y * 4 + 3] = 255;
    }
    ctx.putImageData(this.slice, x, 0);
    ctx.restore();
  }
}

class LogarithmicSpectrogram {
  constructor() {
    this.setRange(1, 1000);
  }

  setRange(low, high) {
    if (low <= 0)
      throw new RangeError('Range low limit must be greater than 0');
    if (low >= high)
      throw new Error('Range low limit less than high limit');
    this.low = low;
    this.high = high;
    // to map [low, high] logarithmically to [0, M]:
    // mapped = log(x - low + 1) / log(high - low + 1) * M
    //
    // to map back:
    // unmapped = base^(x * log(high - low + 1) / M) + low - 1
    // base is the base of used logarithm (10 gor lo10, e for ln, it does not really matter)
    this.invLog = 1 / Math.log(high - low + 1);
  }

  drawSlice(ctx, frequencyData, intensityToColor,) {
    ctx.save();
    ctx.translate(ctx.canvas.width - 0.5, -0.5);
    const height = ctx.canvas.height;
    let y = 0;
    let currentIntensity = 0;
    for (let i = Math.trunc(this.low); i != frequencyData.length; ++i) {
      currentIntensity = Math.max(currentIntensity, frequencyData[i]);
      const coord = Math.log(i - this.low + 1) * this.invLog * height;
      if (coord >= y + 1) {
        console.log();
        const [r, g, b] = intensityToColor(currentIntensity);
        ctx.fillStyle = `rgb(${r},${g},${b})`;
        const c = Math.trunc(coord);
        ctx.fillRect(0, height - c, 1, c - y);
        y = c;
        currentIntensity = 0;
      }
    }
    ctx.restore();
  }
}

const harmonicsCoefficients = [0, 1, 0, 0.15, 0.03, 0.02, 0, 0, 0.01];
const distortionCurve = new Float32Array(1024);
for (let i = 0; i != distortionCurve.length; ++i) {
  const x = i / (distortionCurve.length - 1) * 2 - 1;
  const t = [1, x];
  for (let n = 2; n != harmonicsCoefficients.length; ++n)
    t[n] = 2 * x * t[n - 1] - t[n - 2];
  let norm = 0;
  distortionCurve[i] = t.reduce((sum, value, index) => {
    norm += harmonicsCoefficients[index];
    return sum + value * harmonicsCoefficients[index];
  }, 0);
  distortionCurve[i] /= norm;
}

function scheduleVolumeChange(gainNode, newValue, duration, start = 0) {
  gainNode.gain.cancelScheduledValues(gainNode.context.currentTime);
  const oldValue = gainNode.gain.value;
  start += gainNode.context.currentTime;
  gainNode.gain.linearRampToValueAtTime(oldValue * 0.99 + newValue * 0.01, start + 0.1 * duration);
  gainNode.gain.linearRampToValueAtTime(oldValue * 0.95 + newValue * 0.05, start + 0.2 * duration);
  gainNode.gain.linearRampToValueAtTime(oldValue * 0.85 + newValue * 0.15, start + 0.3 * duration);
  gainNode.gain.linearRampToValueAtTime(oldValue * 0.15 + newValue * 0.85, start + 0.7 * duration);
  gainNode.gain.linearRampToValueAtTime(oldValue * 0.05 + newValue * 0.95, start + 0.8 * duration);
  gainNode.gain.linearRampToValueAtTime(oldValue * 0.01 + newValue * 0.99, start + 0.9 * duration);
  gainNode.gain.linearRampToValueAtTime(newValue, start + duration);
}


(() => {
  const sampleRatePicker = document.getElementById('sampleRate');
  setDefaulSampleRateOption(sampleRatePicker);

  const distortionOption = document.getElementById('distortion');

  const spectrumCtx = document.getElementById('spectrum').getContext('2d', { alpha: false });
  const spectrumMidFrequency = document.getElementById('spectrumMidFrequency');
  const spectrumMaxFrequency = document.getElementById('spectrumMaxFrequency');


  const ctx = document.getElementById('spectrogram').getContext('2d', { alpha: false });
  const spectrogramMidFrequency = document.getElementById('spectrogramMidFrequency');
  const spectrogramMaxFrequency = document.getElementById('spectrogramMaxFrequency');

  const container = document.getElementById('spectrogramContainer');

  function setMidFrequencyLabel(value) {
    const text = formatHz(value);
    spectrumMidFrequency.textContent = text;
    spectrogramMidFrequency.textContent = text;
  }

  function setMaxFrequencyLabel(value) {
    const text = formatHz(value);
    spectrumMaxFrequency.textContent = text;
    spectrogramMaxFrequency.textContent = text;
  }

  let colorScheme;

  function clearCanvases() {
    fillCanvas(spectrumCtx, colorScheme.baseColor);
    fillCanvas(ctx, colorScheme.baseColor);
  }

  const blackOnWhiteColorOption = document.querySelector('#colorScheme .black');
  colorScheme = setupColorSchemeChooser('#colorScheme', spectrumCtx, (scheme) => {
    colorScheme = scheme;
    clearCanvases();
    if (blackOnWhiteColorOption.checked)
      container.classList.add('blackOnWhite');
    else
      container.classList.remove('blackOnWhite');
  });
  clearCanvases();

  function drawSpectrum(ctx, frequencyData, gradient) {
    ctx.save();
    fillCanvas(ctx, colorScheme.baseColor);
    ctx.setTransform(1, 0, 0, -1, 0, ctx.canvas.height);
    const xScale = ctx.canvas.width / frequencyData.length;
    const yScale = ctx.canvas.height / 255;
    ctx.strokeStyle = gradient;
    ctx.fillStyle = gradient;
    ctx.lineWidth = 2;
    ctx.lineCap = 'round';
    ctx.lineJoin = 'round';
    ctx.beginPath();
    ctx.moveTo(0, -1);
    for (let i = 0; i != frequencyData.length; ++i)
      ctx.lineTo(i * xScale, frequencyData[i] * yScale);
    ctx.lineTo(ctx.canvas.width, -1);
    ctx.globalAlpha = 0.5;
    ctx.fill();
    ctx.globalAlpha = 1;
    ctx.stroke();
    ctx.restore();
  }

  const playButton = document.getElementById('playButton');

  function createAudioCtx() {
    const frequencySlider = document.getElementById('frequencySlider');
    function getFrequency() {
      const baseFrequency = 55;
      return baseFrequency * Math.pow(2, Number.parseFloat(frequencySlider.value) / 12);
    };

    const volumeSlider = document.getElementById('volumeSlider');
    const getVolume = () => Number.parseFloat(volumeSlider.value);

    let audioCtx = null;
    let oscillator = null;
    let analyser = null;
    let distortion = null;
    let gain = null;


    function create() {
      const sampleRate = sampleRatePicker.value;
      const options = {};
      if (sampleRate != '')
        options.sampleRate = Number.parseInt(sampleRate);
      audioCtx = new AudioContext(options);

      setMidFrequencyLabel(audioCtx.sampleRate / 4);
      setMaxFrequencyLabel(audioCtx.sampleRate / 2);

      oscillator = audioCtx.createOscillator();
      oscillator.type = 'sine';
      oscillator.frequency.value = getFrequency();
      oscillator.start();

      analyser = audioCtx.createAnalyser();
      analyser.fftSize = 2048;//16384;
      analyser.minDecibels = -90;
      //analyser.maxDecibels = 0;
      analyser.smoothingTimeConstant = 0;

      distortion = audioCtx.createWaveShaper();
      if (distortionOption.checked)
        distortion.curve = distortionCurve;
      distortion.oversample = 'none';

      gain = audioCtx.createGain();
      gain.gain.value = getVolume();

      oscillator.connect(distortion);
      distortion.connect(analyser);
      analyser.connect(gain);
      gain.connect(audioCtx.destination);
    }
    create();

    distortionOption.addEventListener('change', () => {
      distortion.curve = distortionOption.checked ? distortionCurve : null;
    });

    const frequencyDisplay = document.getElementById('frequencyDisplay');
    frequencySlider.addEventListener('input', () => {
      const f = getFrequency();
      if (f <= audioCtx.sampleRate / 2) {
        frequencyDisplay.textContent = f.toFixed() + ' Hz';
        oscillator.frequency.linearRampToValueAtTime(f, audioCtx.currentTime + 0.05);
      }
      else {
        frequencyDisplay.textContent = 'too high for sample rate';
      }
    });

    volumeSlider.addEventListener('input', getThrottledCallback(() => {
      scheduleVolumeChange(gain, getVolume(), 0.05);
    }, 50));

    const spectrogram = new LinearSpectrogram(ctx);

    //const logSpectrogram = new LogarithmicSpectrogram();
    //logSpectrogram.setRange(1, analyser.fftSize / 2 - 1);

    function animateFrame() {
      const shiftWidth = ctx.canvas.width - 1;
      ctx.drawImage(ctx.canvas, 1, 0, shiftWidth, ctx.canvas.height, 0, 0, shiftWidth, ctx.canvas.height);
      const data = new Uint8Array(analyser.frequencyBinCount);
      analyser.getByteFrequencyData(data);
      spectrogram.drawSlice(ctx, ctx.canvas.width - 1, data, colorScheme.intensity);
      //logSpectrogram.drawSlice(ctx, data, colorScheme.intensity);

      drawSpectrum(spectrumCtx, data, colorScheme.gradient);
    }

    const framePeriod = 1000 / 60; // 60 fps
    let id = setInterval(animateFrame, framePeriod)

    playButton.removeEventListener('click', createAudioCtx);

    let suspendID = null;
    playButton.addEventListener('click', () => {
      if (suspendID !== null) {
        clearTimeout(suspendID);
        suspendID = null;
      }
      scheduleVolumeChange(gain, getVolume(), 0.05);
      if (audioCtx.state == 'suspended')
        audioCtx.resume().then(() => {
          id = setInterval(animateFrame, framePeriod)
        });
    });
    document.getElementById('stopButton').addEventListener('click', () => {
      if (suspendID !== null)
        return;

      clearInterval(id);
      scheduleVolumeChange(gain, 0, 0.05);
      suspendID = setTimeout(() => {
        audioCtx.suspend();
      }, 100);
    });

    sampleRatePicker.addEventListener('change', () => {
      const wasPlaying = audioCtx.state == 'running';
      audioCtx.close();
      clearCanvases();
      create();
      if (!wasPlaying)
        audioCtx.suspend();
    });
  }

  playButton.addEventListener('click', createAudioCtx);
})();


(() => {
  const sampleRatePicker = document.getElementById('oversamplingSampleRate');
  const getSampleRate = () => Number.parseInt(sampleRatePicker.value);

  const distortionOption = document.getElementById('oversamplingDistortion');

  const noOversampleOption = document.getElementById('noOversample');
  const oversample2xOption = document.getElementById('2xOversample');
  const oversample4xOption = document.getElementById('4xOversample');
  function getOversampling() {
    if (oversample2xOption.checked)
      return 2;
    else if (oversample4xOption.checked)
      return 4;
    return 1;
  }


  const ctx = document.getElementById('oversamplingSpectrogram').getContext('2d', { alpha: false });
  const spectrogramMidFrequency = document.getElementById('oversamplingSpectrogramMidFrequency');
  const spectrogramMaxFrequency = document.getElementById('oversamplingSpectrogramMaxFrequency');

  const oversampledCtx = document.getElementById('oversampledSpectrogram').getContext('2d', { alpha: false });
  const oversampledSpectrogramMidFrequency = document.getElementById('oversampledSpectrogramMidFrequency');
  const oversampledSpectrogramMaxFrequency = document.getElementById('oversampledSpectrogramMaxFrequency');

  const filterDisplayCtx = document.getElementById('oversamplingFilterDisplay').getContext('2d', { alpha: false });
  const filterDisplayMidFrequency = document.getElementById('oversamplingFilterDisplayMidFrequency');
  const filterDisplayMaxFrequency = document.getElementById('oversamplingFilterDisplayMaxFrequency');

  const container = document.getElementById('oversamplingContainer');

  const frequencySlider = document.getElementById('oversamplingFrequencySlider');

  function setFrecuencyLabels(sampleRate) {
    spectrogramMidFrequency.textContent = formatHz(sampleRate / 4);
    spectrogramMaxFrequency.textContent = formatHz(sampleRate / 2);
  }

  function setOversampledFrequencyLabels(sampleRate, oversample) {
    const oversampledRate = sampleRate * oversample;
    const oversampledMid = formatHz(oversampledRate / 4);
    const oversampledMax = formatHz(oversampledRate / 2)
    oversampledSpectrogramMidFrequency.textContent = oversampledMid;
    oversampledSpectrogramMaxFrequency.textContent = oversampledMax;
    filterDisplayMidFrequency.textContent = oversampledMid;
    filterDisplayMaxFrequency.textContent = oversampledMax;
  };

  let colorScheme;

  function clearCanvases() {
    fillCanvas(ctx, colorScheme.baseColor);
    fillCanvas(oversampledCtx, colorScheme.baseColor);
  }

  let render = () => {
    ctx.save();
    ctx.fillStyle = 'gray';
    ctx.textBaseline = 'top';
    ctx.textAlign = 'right';
    ctx.font = '16px sans-serif';
    ctx.fillText('loading module...', ctx.canvas.width - 5, 5);
    ctx.restore();
  };

  const drawFrequencyResponse = () => {
    fillCanvas(filterDisplayCtx, colorScheme.baseColor);
    const f0 = Number.parseFloat(frequencySlider.value);
    const w0 = 2 * Math.PI * f0;
    const cosw0 = Math.cos(w0);
    const sinw0 = Math.sin(w0);
    function getGain(Q) {
      const alpha = sinw0 / (2 * Q);
      const invA0 = 1 / (1 + alpha);  // norm to a0 == 1
      const a1 = -2 * cosw0 * invA0;
      const a2 = (1 - alpha) * invA0;
      const b0 = (1 - cosw0) / 2 * invA0;
      const b1 = (1 - cosw0) * invA0;
      const b2 = (1 - cosw0) / 2 * invA0;

      const sqr = x => x * x;

      // numerator = 16*b0*b2*phi^2 - 4*(b0*b1 + 4*b0*b2 + b1*b2)*phi + (b0 + b1 + b2)^2
      // denominator = 16*a2*phi^2 - 4*(a1 + 4*a2 + a1*a2)*phi + (a1 + a2)^2
      // f(phi)=(numerator/denimunator)^0.5
      const c2 = 16 * b0 * b2;
      const c1 = -4 * (b0 * b1 + 4 * b0 * b2 + b1 * b2);
      const c0 = sqr(b0 + b1 + b2);
      const d2 = 16 * a2;
      const d1 = -4 * (a1 + 4 * a2 + a1 * a2);
      const d0 = sqr(1 + a1 + a2);
      return function (w) {
        const phi = sqr(Math.sin(Math.PI * w * 0.5));
        return 10 * (Math.log10((c2 * phi + c1) * phi + c0) - Math.log10((d2 * phi + d1) * phi + d0));
      };
    }

    const f1 = getGain(0.54119610);
    const f2 = getGain(1.3065630);
    const f = w => f1(w) + f2(w);

    const dBMin = -90;
    const dBMax = 0;
    const deltaDb = 1 / (dBMax - dBMin);
    const dBMap = dB => deltaDb * (dB - dBMin);

    filterDisplayCtx.save();
    filterDisplayCtx.setTransform(1, 0, 0, -1, 0, filterDisplayCtx.canvas.height + 10);
    filterDisplayCtx.strokeStyle = colorScheme.gradient;
    filterDisplayCtx.fillStyle = colorScheme.gradient;
    filterDisplayCtx.lineWidth = 2;
    filterDisplayCtx.lineCap = 'round';
    filterDisplayCtx.lineJoin = 'round';
    filterDisplayCtx.beginPath();
    filterDisplayCtx.moveTo(0, -1);
    const pointCount = 300;
    const pixelsPerPoint = filterDisplayCtx.canvas.width / (pointCount - 1);
    for (let i = 0; i != pointCount; ++i) {
      const dB = f(i / pointCount);
      filterDisplayCtx.lineTo(i * pixelsPerPoint, Math.max(dBMap(dB) * filterDisplayCtx.canvas.height, -1));
    }
    filterDisplayCtx.lineTo(filterDisplayCtx.canvas.width, -1);
    filterDisplayCtx.globalAlpha = 0.3;
    filterDisplayCtx.fill();
    filterDisplayCtx.globalAlpha = 1;
    filterDisplayCtx.stroke();

    filterDisplayCtx.strokeStyle = '#7777';
    filterDisplayCtx.lineWidth = 1;
    filterDisplayCtx.beginPath();
    const x0 = Math.round(2 * f0 * filterDisplayCtx.canvas.width) + 0.5;
    filterDisplayCtx.moveTo(x0, 0);
    filterDisplayCtx.lineTo(x0, dBMap(-3) * filterDisplayCtx.canvas.height);
    filterDisplayCtx.stroke();
    filterDisplayCtx.restore();
  };

  const blackOnWhiteColorOption = document.querySelector('#oversamplingColorScheme .black');
  colorScheme = setupColorSchemeChooser('#oversamplingColorScheme', filterDisplayCtx, (scheme) => {
    colorScheme = scheme;
    clearCanvases();
    render();
    drawFrequencyResponse();
    if (blackOnWhiteColorOption.checked)
      container.classList.add('blackOnWhite');
    else
      container.classList.remove('blackOnWhite');
  });
  clearCanvases();
  render();
  drawFrequencyResponse();

  Module.postRun.push(() => {
    let oversampling = new Module.Oversampling();
    oversampling.minDecibels = -90;

    const sampleRate = getSampleRate();
    const oversample = getOversampling();

    setFrecuencyLabels(sampleRate);
    setOversampledFrequencyLabels(sampleRate, oversample);
    oversampling.sampleRate = sampleRate;

    oversampling.oversample = oversample;
    oversampling.fftSize = 2048;

    const oversamplingSpectrogramRenderingLabel = document.getElementById('oversamplingSpectrogramRenderingLabel');
    const oversampledSpectrogramRenderingLabel = document.getElementById('oversampledSpectrogramRenderingLabel');
    const playButton = document.getElementById('oversamplingPlayButton');

    let bufferIsDirty = true;

    const spectrogram = new LinearSpectrogram(ctx);
    const oversampledSpectrogram = new LinearSpectrogram(oversampledCtx);

    fillCanvas(ctx, colorScheme.baseColor);

    let renderStep;
    let renderStartTime;
    render = () => {
      if (debug) {
        renderStartTime = performance.now();
        console.log(`restarting with
Fs: ${oversampling.sampleRate} Hz
sweep start: ${oversampling.startFrequency} Hz
      end: ${Math.min(oversampling.endFrequency, oversampling.sampleRate / 2)} (set ${oversampling.endFrequency}) Hz
fft size: ${oversampling.fftSize}
oversample: ${oversampling.oversample}
cutoff frequency: ${oversampling.cutoffFrequency} Hz
harmonics coefficients: ${oversampling.harmonicsCoefficients}`);
      }
      let timestamp = performance.now();
      oversamplingSpectrogramRenderingLabel.classList.remove('hidden');
      oversampledSpectrogramRenderingLabel.classList.remove('hidden');
      playButton.setAttribute('disabled', null);
      renderStep = 0;
      oversampling.restartAndRenderStep();
      spectrogram.drawSlice(ctx, 0, oversampling.stepFrequencyData, colorScheme.intensity);
      oversampledSpectrogram.drawSlice(oversampledCtx, 0, oversampling.stepOversampledFrequencyData, colorScheme.intensity);
      function renderChunk() {
        const renderChunkSize = 25;
        for (let i = 0; i != renderChunkSize; ++i) {
          if (++renderStep >= ctx.canvas.width) {
            bufferIsDirty = true;
            oversamplingSpectrogramRenderingLabel.classList.add('hidden');
            oversampledSpectrogramRenderingLabel.classList.add('hidden');
            playButton.removeAttribute('disabled');
            if (debug) {
              console.log(`Rendered in ${performance.now() - renderStartTime} ms`);
            }
            return;
          }
          oversampling.renderNextStep();
          spectrogram.drawSlice(ctx, renderStep, oversampling.stepFrequencyData, colorScheme.intensity);
          oversampledSpectrogram.drawSlice(oversampledCtx, renderStep, oversampling.stepOversampledFrequencyData, colorScheme.intensity);
        }
        setTimeout(renderChunk);
      }
      renderChunk();
    };

    if (distortionOption.checked)
      oversampling.harmonicsCoefficients = harmonicsCoefficients;
    function resetCutoffFrequency(sampleRate, oversample) {
      oversampling.cutoffFrequency = sampleRate * oversample * Number.parseFloat(frequencySlider.value);
    }

    const volumeSlider = document.getElementById('oversamplingVolumeSlider');
    const getVolume = () => Number.parseFloat(volumeSlider.value);

    let audioCtx = null;
    let audioBuffer = null;
    let gain;

    function createAudioCtx() {
      audioCtx = new AudioContext({ sampleRate: getSampleRate() });
      audioCtx.suspend();
      gain = audioCtx.createGain();
      gain.gain.value = getVolume();
      gain.connect(audioCtx.destination);

      audioBuffer = audioCtx.createBuffer(1, oversampling.audioBufferSampleCount, audioCtx.sampleRate);
    }

    volumeSlider.addEventListener('input', getThrottledCallback(() => {
      scheduleVolumeChange(gain, getVolume(), 0.05);
    }, 50));

    playButton.addEventListener('click', () => {
      if (audioCtx === null)
        createAudioCtx();

      if (audioCtx.state != 'suspended')
        return;

      if (bufferIsDirty) {
        audioBuffer.copyToChannel(oversampling.audioBufferData, 0);
        bufferIsDirty = false;
      }

      gain.gain.cancelScheduledValues(audioCtx.currentTime);
      gain.gain.value = getVolume();

      const audioBufferNode = audioCtx.createBufferSource();
      audioBufferNode.buffer = audioBuffer;
      audioBufferNode.connect(gain);
      audioBufferNode.start();
      audioCtx.resume().then(() => {
        scheduleVolumeChange(gain, 0, 0.05, 1.950);
        setTimeout(() => {
          audioCtx.suspend();
        }, 2050);
      });
    });

    distortionOption.addEventListener('change', () => {
      oversampling.harmonicsCoefficients = distortionOption.checked ? harmonicsCoefficients : [];
      render();
    });

    sampleRatePicker.addEventListener('change', () => {
      const sampleRate = getSampleRate();
      const oversample = getOversampling();
      oversampling.sampleRate = sampleRate;
      resetCutoffFrequency(sampleRate, oversample);
      setFrecuencyLabels(sampleRate);
      setOversampledFrequencyLabels(sampleRate, oversample);
      render();
      if (audioCtx !== null) {
        audioCtx.close();
        audioCtx = null;
      }
    });

    const handleOversamplingChange = () => {
      const sampleRate = getSampleRate();
      const oversample = getOversampling();
      oversampling.oversample = oversample;
      resetCutoffFrequency(sampleRate, oversample);
      setOversampledFrequencyLabels(sampleRate, oversample);
      render();
    };
    noOversampleOption.addEventListener('change', handleOversamplingChange);
    oversample2xOption.addEventListener('change', handleOversamplingChange);
    oversample4xOption.addEventListener('change', handleOversamplingChange);

    const renderThrottled = getThrottledCallback(() => {
      render();
    }, 200, 200);
    oversamplingFrequencySlider.addEventListener('input', getDrawOnNextAnimationFrameThrottled(() => {
      const sampleRate = getSampleRate();
      const oversample = getOversampling();
      resetCutoffFrequency(sampleRate, oversample);
      drawFrequencyResponse();
      renderThrottled();
    }));

    render();
  });
})();