
var f;f||(f=typeof Module !== 'undefined' ? Module : {});var ba={},q;for(q in f)f.hasOwnProperty(q)&&(ba[q]=f[q]);var ca=!1,y=!1,da=!1,ea=!1;ca="object"===typeof window;y="function"===typeof importScripts;da="object"===typeof process&&"object"===typeof process.versions&&"string"===typeof process.versions.node;ea=!ca&&!da&&!y;var z="",fa,ha,ia,ja;
if(da)z=y?require("path").dirname(z)+"/":__dirname+"/",fa=function(a,b){ia||(ia=require("fs"));ja||(ja=require("path"));a=ja.normalize(a);return ia.readFileSync(a,b?null:"utf8")},ha=function(a){a=fa(a,!0);a.buffer||(a=new Uint8Array(a));assert(a.buffer);return a},1<process.argv.length&&process.argv[1].replace(/\\/g,"/"),process.argv.slice(2),"undefined"!==typeof module&&(module.exports=f),process.on("uncaughtException",function(a){throw a;}),process.on("unhandledRejection",A),f.inspect=function(){return"[Emscripten Module object]"};
else if(ea)"undefined"!=typeof read&&(fa=function(a){return read(a)}),ha=function(a){if("function"===typeof readbuffer)return new Uint8Array(readbuffer(a));a=read(a,"binary");assert("object"===typeof a);return a},"undefined"!==typeof print&&("undefined"===typeof console&&(console={}),console.log=print,console.warn=console.error="undefined"!==typeof printErr?printErr:print);else if(ca||y)y?z=self.location.href:"undefined"!==typeof document&&document.currentScript&&(z=document.currentScript.src),z=
0!==z.indexOf("blob:")?z.substr(0,z.lastIndexOf("/")+1):"",fa=function(a){var b=new XMLHttpRequest;b.open("GET",a,!1);b.send(null);return b.responseText},y&&(ha=function(a){var b=new XMLHttpRequest;b.open("GET",a,!1);b.responseType="arraybuffer";b.send(null);return new Uint8Array(b.response)});var ka=f.print||console.log.bind(console),B=f.printErr||console.warn.bind(console);for(q in ba)ba.hasOwnProperty(q)&&(f[q]=ba[q]);ba=null;var la;f.wasmBinary&&(la=f.wasmBinary);var noExitRuntime;
f.noExitRuntime&&(noExitRuntime=f.noExitRuntime);"object"!==typeof WebAssembly&&A("no native wasm support detected");var ma,na=!1;function assert(a,b){a||A("Assertion failed: "+b)}var oa="undefined"!==typeof TextDecoder?new TextDecoder("utf8"):void 0;
function pa(a,b,c){var d=b+c;for(c=b;a[c]&&!(c>=d);)++c;if(16<c-b&&a.subarray&&oa)return oa.decode(a.subarray(b,c));for(d="";b<c;){var e=a[b++];if(e&128){var g=a[b++]&63;if(192==(e&224))d+=String.fromCharCode((e&31)<<6|g);else{var k=a[b++]&63;e=224==(e&240)?(e&15)<<12|g<<6|k:(e&7)<<18|g<<12|k<<6|a[b++]&63;65536>e?d+=String.fromCharCode(e):(e-=65536,d+=String.fromCharCode(55296|e>>10,56320|e&1023))}}else d+=String.fromCharCode(e)}return d}
function qa(a,b,c){var d=D;if(0<c){c=b+c-1;for(var e=0;e<a.length;++e){var g=a.charCodeAt(e);if(55296<=g&&57343>=g){var k=a.charCodeAt(++e);g=65536+((g&1023)<<10)|k&1023}if(127>=g){if(b>=c)break;d[b++]=g}else{if(2047>=g){if(b+1>=c)break;d[b++]=192|g>>6}else{if(65535>=g){if(b+2>=c)break;d[b++]=224|g>>12}else{if(b+3>=c)break;d[b++]=240|g>>18;d[b++]=128|g>>12&63}d[b++]=128|g>>6&63}d[b++]=128|g&63}}d[b]=0}}var ra="undefined"!==typeof TextDecoder?new TextDecoder("utf-16le"):void 0;
function sa(a,b){var c=a>>1;for(var d=c+b/2;!(c>=d)&&ta[c];)++c;c<<=1;if(32<c-a&&ra)return ra.decode(D.subarray(a,c));c="";for(d=0;!(d>=b/2);++d){var e=E[a+2*d>>1];if(0==e)break;c+=String.fromCharCode(e)}return c}function ua(a,b,c){void 0===c&&(c=2147483647);if(2>c)return 0;c-=2;var d=b;c=c<2*a.length?c/2:a.length;for(var e=0;e<c;++e)E[b>>1]=a.charCodeAt(e),b+=2;E[b>>1]=0;return b-d}function va(a){return 2*a.length}
function wa(a,b){for(var c=0,d="";!(c>=b/4);){var e=F[a+4*c>>2];if(0==e)break;++c;65536<=e?(e-=65536,d+=String.fromCharCode(55296|e>>10,56320|e&1023)):d+=String.fromCharCode(e)}return d}function xa(a,b,c){void 0===c&&(c=2147483647);if(4>c)return 0;var d=b;c=d+c-4;for(var e=0;e<a.length;++e){var g=a.charCodeAt(e);if(55296<=g&&57343>=g){var k=a.charCodeAt(++e);g=65536+((g&1023)<<10)|k&1023}F[b>>2]=g;b+=4;if(b+4>c)break}F[b>>2]=0;return b-d}
function ya(a){for(var b=0,c=0;c<a.length;++c){var d=a.charCodeAt(c);55296<=d&&57343>=d&&++c;b+=4}return b}var za,Aa,D,E,ta,F,G,Ba,Ca;function Da(){var a=ma.buffer;za=a;f.HEAP8=Aa=new Int8Array(a);f.HEAP16=E=new Int16Array(a);f.HEAP32=F=new Int32Array(a);f.HEAPU8=D=new Uint8Array(a);f.HEAPU16=ta=new Uint16Array(a);f.HEAPU32=G=new Uint32Array(a);f.HEAPF32=Ba=new Float32Array(a);f.HEAPF64=Ca=new Float64Array(a)}var Ea,Fa=[],Ga=[],Ha=[],Ia=[];function Ja(){var a=f.preRun.shift();Fa.unshift(a)}
var H=0,Ka=null,La=null;f.preloadedImages={};f.preloadedAudios={};function A(a){if(f.onAbort)f.onAbort(a);B(a);na=!0;throw new WebAssembly.RuntimeError("abort("+a+"). Build with -s ASSERTIONS=1 for more info.");}function Ma(a){var b=J;return String.prototype.startsWith?b.startsWith(a):0===b.indexOf(a)}function Na(){return Ma("data:application/octet-stream;base64,")}var J="fft.wasm";if(!Na()){var Oa=J;J=f.locateFile?f.locateFile(Oa,z):z+Oa}
function Pa(){var a=J;try{if(a==J&&la)return new Uint8Array(la);if(ha)return ha(a);throw"both async and sync fetching of the wasm failed";}catch(b){A(b)}}function Qa(){return la||!ca&&!y||"function"!==typeof fetch||Ma("file://")?Promise.resolve().then(function(){return Pa()}):fetch(J,{credentials:"same-origin"}).then(function(a){if(!a.ok)throw"failed to load wasm binary file at '"+J+"'";return a.arrayBuffer()}).catch(function(){return Pa()})}
function Ra(a){for(;0<a.length;){var b=a.shift();if("function"==typeof b)b(f);else{var c=b.pa;"number"===typeof c?void 0===b.ia?Ea.get(c)():Ea.get(c)(b.ia):c(void 0===b.ia?null:b.ia)}}}function Sa(a){this.R=a-16;this.Da=function(b){F[this.R+8>>2]=b};this.Aa=function(b){F[this.R+0>>2]=b};this.Ba=function(){F[this.R+4>>2]=0};this.za=function(){Aa[this.R+12>>0]=0};this.Ca=function(){Aa[this.R+13>>0]=0};this.ta=function(b,c){this.Da(b);this.Aa(c);this.Ba();this.za();this.Ca()}}var Ta=0;
function Ua(a){switch(a){case 1:return 0;case 2:return 1;case 4:return 2;case 8:return 3;default:throw new TypeError("Unknown type size: "+a);}}var Va=void 0;function K(a){for(var b="";D[a];)b+=Va[D[a++]];return b}var L={},M={},Wa={};function Xa(a){if(void 0===a)return"_unknown";a=a.replace(/[^a-zA-Z0-9_]/g,"$");var b=a.charCodeAt(0);return 48<=b&&57>=b?"_"+a:a}
function Ya(a,b){a=Xa(a);return(new Function("body","return function "+a+'() {\n    "use strict";    return body.apply(this, arguments);\n};\n'))(b)}function Za(a){var b=Error,c=Ya(a,function(d){this.name=a;this.message=d;d=Error(d).stack;void 0!==d&&(this.stack=this.toString()+"\n"+d.replace(/^Error(:[^\n]*)?\n/,""))});c.prototype=Object.create(b.prototype);c.prototype.constructor=c;c.prototype.toString=function(){return void 0===this.message?this.name:this.name+": "+this.message};return c}
var N=void 0;function O(a){throw new N(a);}var $a=void 0;function ab(a){throw new $a(a);}function P(a,b,c){function d(h){h=c(h);h.length!==a.length&&ab("Mismatched type converter count");for(var l=0;l<a.length;++l)Q(a[l],h[l])}a.forEach(function(h){Wa[h]=b});var e=Array(b.length),g=[],k=0;b.forEach(function(h,l){M.hasOwnProperty(h)?e[l]=M[h]:(g.push(h),L.hasOwnProperty(h)||(L[h]=[]),L[h].push(function(){e[l]=M[h];++k;k===g.length&&d(e)}))});0===g.length&&d(e)}
function Q(a,b,c){c=c||{};if(!("argPackAdvance"in b))throw new TypeError("registerType registeredInstance requires argPackAdvance");var d=b.name;a||O('type "'+d+'" must have a positive integer typeid pointer');if(M.hasOwnProperty(a)){if(c.sa)return;O("Cannot register type '"+d+"' twice")}M[a]=b;delete Wa[a];L.hasOwnProperty(a)&&(b=L[a],delete L[a],b.forEach(function(e){e()}))}function bb(a){return{count:a.count,ba:a.ba,da:a.da,R:a.R,T:a.T,U:a.U,V:a.V}}
function cb(a){O(a.P.T.S.name+" instance already deleted")}var db=!1;function eb(){}function fb(a){--a.count.value;0===a.count.value&&(a.U?a.V.aa(a.U):a.T.S.aa(a.R))}function gb(a){if("undefined"===typeof FinalizationGroup)return gb=function(b){return b},a;db=new FinalizationGroup(function(b){for(var c=b.next();!c.done;c=b.next())c=c.value,c.R?fb(c):console.warn("object already deleted: "+c.R)});gb=function(b){db.register(b,b.P,b.P);return b};eb=function(b){db.unregister(b.P)};return gb(a)}
var hb=void 0,ib=[];function jb(){for(;ib.length;){var a=ib.pop();a.P.ba=!1;a["delete"]()}}function R(){}var kb={};function lb(a,b,c){if(void 0===a[b].X){var d=a[b];a[b]=function(){a[b].X.hasOwnProperty(arguments.length)||O("Function '"+c+"' called with an invalid number of arguments ("+arguments.length+") - expects one of ("+a[b].X+")!");return a[b].X[arguments.length].apply(this,arguments)};a[b].X=[];a[b].X[d.fa]=d}}
function mb(a,b){f.hasOwnProperty(a)?(O("Cannot register public name '"+a+"' twice"),lb(f,a,a),f.hasOwnProperty(void 0)&&O("Cannot register multiple overloads of a function with the same number of arguments (undefined)!"),f[a].X[void 0]=b):f[a]=b}function nb(a,b,c,d,e,g,k,h){this.name=a;this.constructor=b;this.Z=c;this.aa=d;this.W=e;this.qa=g;this.ea=k;this.oa=h;this.wa=[]}
function ob(a,b,c){for(;b!==c;)b.ea||O("Expected null or instance of "+c.name+", got an instance of "+b.name),a=b.ea(a),b=b.W;return a}function pb(a,b){if(null===b)return this.ja&&O("null is not a valid "+this.name),0;b.P||O('Cannot pass "'+U(b)+'" as a '+this.name);b.P.R||O("Cannot pass deleted object as a pointer of type "+this.name);return ob(b.P.R,b.P.T.S,this.S)}
function qb(a,b){if(null===b){this.ja&&O("null is not a valid "+this.name);if(this.ha){var c=this.xa();null!==a&&a.push(this.aa,c);return c}return 0}b.P||O('Cannot pass "'+U(b)+'" as a '+this.name);b.P.R||O("Cannot pass deleted object as a pointer of type "+this.name);!this.ga&&b.P.T.ga&&O("Cannot convert argument of type "+(b.P.V?b.P.V.name:b.P.T.name)+" to parameter type "+this.name);c=ob(b.P.R,b.P.T.S,this.S);if(this.ha)switch(void 0===b.P.U&&O("Passing raw pointer to smart pointer is illegal"),
this.Ea){case 0:b.P.V===this?c=b.P.U:O("Cannot convert argument of type "+(b.P.V?b.P.V.name:b.P.T.name)+" to parameter type "+this.name);break;case 1:c=b.P.U;break;case 2:if(b.P.V===this)c=b.P.U;else{var d=b.clone();c=this.ya(c,V(function(){d["delete"]()}));null!==a&&a.push(this.aa,c)}break;default:O("Unsupporting sharing policy")}return c}
function rb(a,b){if(null===b)return this.ja&&O("null is not a valid "+this.name),0;b.P||O('Cannot pass "'+U(b)+'" as a '+this.name);b.P.R||O("Cannot pass deleted object as a pointer of type "+this.name);b.P.T.ga&&O("Cannot convert argument of type "+b.P.T.name+" to parameter type "+this.name);return ob(b.P.R,b.P.T.S,this.S)}function sb(a){return this.fromWireType(G[a>>2])}function tb(a,b,c){if(b===c)return a;if(void 0===c.W)return null;a=tb(a,b,c.W);return null===a?null:c.oa(a)}var ub={};
function vb(a,b){for(void 0===b&&O("ptr should not be undefined");a.W;)b=a.ea(b),a=a.W;return ub[b]}function wb(a,b){b.T&&b.R||ab("makeClassHandle requires ptr and ptrType");!!b.V!==!!b.U&&ab("Both smartPtrType and smartPtr must be specified");b.count={value:1};return gb(Object.create(a,{P:{value:b}}))}
function W(a,b,c,d){this.name=a;this.S=b;this.ja=c;this.ga=d;this.ha=!1;this.aa=this.ya=this.xa=this.ma=this.Ea=this.va=void 0;void 0!==b.W?this.toWireType=qb:(this.toWireType=d?pb:rb,this.Y=null)}function zb(a,b){f.hasOwnProperty(a)||ab("Replacing nonexistant public symbol");f[a]=b;f[a].fa=void 0}
function Ab(a,b){assert(0<=a.indexOf("j"),"getDynCaller should only be called with i64 sigs");var c=[];return function(){c.length=arguments.length;for(var d=0;d<arguments.length;d++)c[d]=arguments[d];var e;-1!=a.indexOf("j")?e=c&&c.length?f["dynCall_"+a].apply(null,[b].concat(c)):f["dynCall_"+a].call(null,b):e=Ea.get(b).apply(null,c);return e}}function X(a,b){a=K(a);var c=-1!=a.indexOf("j")?Ab(a,b):Ea.get(b);"function"!==typeof c&&O("unknown function pointer with signature "+a+": "+b);return c}
var Bb=void 0;function Cb(a){a=Db(a);var b=K(a);Y(a);return b}function Eb(a,b){function c(g){e[g]||M[g]||(Wa[g]?Wa[g].forEach(c):(d.push(g),e[g]=!0))}var d=[],e={};b.forEach(c);throw new Bb(a+": "+d.map(Cb).join([", "]));}function Fb(a,b){for(var c=[],d=0;d<a;d++)c.push(F[(b>>2)+d]);return c}function Gb(a){for(;a.length;){var b=a.pop();a.pop()(b)}}
function Hb(a){var b=Function;if(!(b instanceof Function))throw new TypeError("new_ called with constructor type "+typeof b+" which is not a function");var c=Ya(b.name||"unknownFunctionName",function(){});c.prototype=b.prototype;c=new c;a=b.apply(c,a);return a instanceof Object?a:c}
function Ib(a,b,c){a instanceof Object||O(c+' with invalid "this": '+a);a instanceof b.S.constructor||O(c+' incompatible with "this" of type '+a.constructor.name);a.P.R||O("cannot call emscripten binding method "+c+" on deleted object");return ob(a.P.R,a.P.T.S,b.S)}var Jb=[],Z=[{},{value:void 0},{value:null},{value:!0},{value:!1}];function Kb(a){4<a&&0===--Z[a].ka&&(Z[a]=void 0,Jb.push(a))}
function V(a){switch(a){case void 0:return 1;case null:return 2;case !0:return 3;case !1:return 4;default:var b=Jb.length?Jb.pop():Z.length;Z[b]={ka:1,value:a};return b}}function U(a){if(null===a)return"null";var b=typeof a;return"object"===b||"array"===b||"function"===b?a.toString():""+a}function Lb(a,b){switch(b){case 2:return function(c){return this.fromWireType(Ba[c>>2])};case 3:return function(c){return this.fromWireType(Ca[c>>3])};default:throw new TypeError("Unknown float type: "+a);}}
function Mb(a,b,c){switch(b){case 0:return c?function(d){return Aa[d]}:function(d){return D[d]};case 1:return c?function(d){return E[d>>1]}:function(d){return ta[d>>1]};case 2:return c?function(d){return F[d>>2]}:function(d){return G[d>>2]};default:throw new TypeError("Unknown integer type: "+a);}}function Nb(a){a||O("Cannot use deleted val. handle = "+a);return Z[a].value}function Ob(a,b){var c=M[a];void 0===c&&O(b+" has unknown type "+Cb(a));return c}var Pb={};
function Qb(a){var b=Pb[a];return void 0===b?K(a):b}var Rb=[];function Sb(a){var b=Rb.length;Rb.push(a);return b}function Tb(a,b){for(var c=Array(a),d=0;d<a;++d)c[d]=Ob(F[(b>>2)+d],"parameter "+d);return c}for(var Ub=[null,[],[]],Vb=Array(256),Wb=0;256>Wb;++Wb)Vb[Wb]=String.fromCharCode(Wb);Va=Vb;N=f.BindingError=Za("BindingError");$a=f.InternalError=Za("InternalError");
R.prototype.isAliasOf=function(a){if(!(this instanceof R&&a instanceof R))return!1;var b=this.P.T.S,c=this.P.R,d=a.P.T.S;for(a=a.P.R;b.W;)c=b.ea(c),b=b.W;for(;d.W;)a=d.ea(a),d=d.W;return b===d&&c===a};R.prototype.clone=function(){this.P.R||cb(this);if(this.P.da)return this.P.count.value+=1,this;var a=gb(Object.create(Object.getPrototypeOf(this),{P:{value:bb(this.P)}}));a.P.count.value+=1;a.P.ba=!1;return a};
R.prototype["delete"]=function(){this.P.R||cb(this);this.P.ba&&!this.P.da&&O("Object already scheduled for deletion");eb(this);fb(this.P);this.P.da||(this.P.U=void 0,this.P.R=void 0)};R.prototype.isDeleted=function(){return!this.P.R};R.prototype.deleteLater=function(){this.P.R||cb(this);this.P.ba&&!this.P.da&&O("Object already scheduled for deletion");ib.push(this);1===ib.length&&hb&&hb(jb);this.P.ba=!0;return this};W.prototype.ra=function(a){this.ma&&(a=this.ma(a));return a};
W.prototype.la=function(a){this.aa&&this.aa(a)};W.prototype.argPackAdvance=8;W.prototype.readValueFromPointer=sb;W.prototype.deleteObject=function(a){if(null!==a)a["delete"]()};
W.prototype.fromWireType=function(a){function b(){return this.ha?wb(this.S.Z,{T:this.va,R:c,V:this,U:a}):wb(this.S.Z,{T:this,R:a})}var c=this.ra(a);if(!c)return this.la(a),null;var d=vb(this.S,c);if(void 0!==d){if(0===d.P.count.value)return d.P.R=c,d.P.U=a,d.clone();d=d.clone();this.la(a);return d}d=this.S.qa(c);d=kb[d];if(!d)return b.call(this);d=this.ga?d.na:d.pointerType;var e=tb(c,this.S,d.S);return null===e?b.call(this):this.ha?wb(d.S.Z,{T:d,R:e,V:this,U:a}):wb(d.S.Z,{T:d,R:e})};
f.getInheritedInstanceCount=function(){return Object.keys(ub).length};f.getLiveInheritedInstances=function(){var a=[],b;for(b in ub)ub.hasOwnProperty(b)&&a.push(ub[b]);return a};f.flushPendingDeletes=jb;f.setDelayFunction=function(a){hb=a;ib.length&&hb&&hb(jb)};Bb=f.UnboundTypeError=Za("UnboundTypeError");f.count_emval_handles=function(){for(var a=0,b=5;b<Z.length;++b)void 0!==Z[b]&&++a;return a};f.get_first_emval=function(){for(var a=5;a<Z.length;++a)if(void 0!==Z[a])return Z[a];return null};Ga.push({pa:function(){Xb()}});
var Zb={q:function(a){return Yb(a+16)+16},F:function(a,b,c){(new Sa(a)).ta(b,c);Ta++;throw a;},z:function(a,b,c,d,e){var g=Ua(c);b=K(b);Q(a,{name:b,fromWireType:function(k){return!!k},toWireType:function(k,h){return h?d:e},argPackAdvance:8,readValueFromPointer:function(k){if(1===c)var h=Aa;else if(2===c)h=E;else if(4===c)h=F;else throw new TypeError("Unknown boolean type size: "+b);return this.fromWireType(h[k>>g])},Y:null})},j:function(a,b,c,d,e,g,k,h,l,m,n,r,u){n=K(n);g=X(e,g);h&&(h=X(k,h));m&&
(m=X(l,m));u=X(r,u);var x=Xa(n);mb(x,function(){Eb("Cannot construct "+n+" due to unbound types",[d])});P([a,b,c],d?[d]:[],function(t){t=t[0];if(d){var v=t.S;var p=v.Z}else p=R.prototype;t=Ya(x,function(){if(Object.getPrototypeOf(this)!==C)throw new N("Use 'new' to construct "+n);if(void 0===w.$)throw new N(n+" has no accessible constructor");var S=w.$[arguments.length];if(void 0===S)throw new N("Tried to invoke ctor of "+n+" with invalid number of parameters ("+arguments.length+") - expected ("+
Object.keys(w.$).toString()+") parameters instead!");return S.apply(this,arguments)});var C=Object.create(p,{constructor:{value:t}});t.prototype=C;var w=new nb(n,t,C,u,v,g,h,m);v=new W(n,w,!0,!1);p=new W(n+"*",w,!1,!1);var I=new W(n+" const*",w,!1,!0);kb[a]={pointerType:p,na:I};zb(x,t);return[v,p,I]})},i:function(a,b,c,d,e,g){assert(0<b);var k=Fb(b,c);e=X(d,e);var h=[g],l=[];P([],[a],function(m){m=m[0];var n="constructor "+m.name;void 0===m.S.$&&(m.S.$=[]);if(void 0!==m.S.$[b-1])throw new N("Cannot register multiple constructors with identical number of parameters ("+
(b-1)+") for class '"+m.name+"'! Overload resolution is currently only performed using the parameter count, not actual type info!");m.S.$[b-1]=function(){Eb("Cannot construct "+m.name+" due to unbound types",k)};P([],k,function(r){m.S.$[b-1]=function(){arguments.length!==b-1&&O(n+" called with "+arguments.length+" arguments, expected "+(b-1));l.length=0;h.length=b;for(var u=1;u<b;++u)h[u]=r[u].toWireType(l,arguments[u-1]);u=e.apply(null,h);Gb(l);return r[0].fromWireType(u)};return[]});return[]})},
e:function(a,b,c,d,e,g,k,h){var l=Fb(c,d);b=K(b);g=X(e,g);P([],[a],function(m){function n(){Eb("Cannot call "+r+" due to unbound types",l)}m=m[0];var r=m.name+"."+b;h&&m.S.wa.push(b);var u=m.S.Z,x=u[b];void 0===x||void 0===x.X&&x.className!==m.name&&x.fa===c-2?(n.fa=c-2,n.className=m.name,u[b]=n):(lb(u,b,r),u[b].X[c-2]=n);P([],l,function(t){var v=r,p=m,C=g,w=t.length;2>w&&O("argTypes array size mismatch! Must at least get return value and 'this' types!");var I=null!==t[1]&&null!==p,S=!1;for(p=1;p<
t.length;++p)if(null!==t[p]&&void 0===t[p].Y){S=!0;break}var xb="void"!==t[0].name,T="",aa="";for(p=0;p<w-2;++p)T+=(0!==p?", ":"")+"arg"+p,aa+=(0!==p?", ":"")+"arg"+p+"Wired";v="return function "+Xa(v)+"("+T+") {\nif (arguments.length !== "+(w-2)+") {\nthrowBindingError('function "+v+" called with ' + arguments.length + ' arguments, expected "+(w-2)+" args!');\n}\n";S&&(v+="var destructors = [];\n");var yb=S?"destructors":"null";T="throwBindingError invoker fn runDestructors retType classParam".split(" ");
C=[O,C,k,Gb,t[0],t[1]];I&&(v+="var thisWired = classParam.toWireType("+yb+", this);\n");for(p=0;p<w-2;++p)v+="var arg"+p+"Wired = argType"+p+".toWireType("+yb+", arg"+p+"); // "+t[p+2].name+"\n",T.push("argType"+p),C.push(t[p+2]);I&&(aa="thisWired"+(0<aa.length?", ":"")+aa);v+=(xb?"var rv = ":"")+"invoker(fn"+(0<aa.length?", ":"")+aa+");\n";if(S)v+="runDestructors(destructors);\n";else for(p=I?1:2;p<t.length;++p)w=1===p?"thisWired":"arg"+(p-2)+"Wired",null!==t[p].Y&&(v+=w+"_dtor("+w+"); // "+t[p].name+
"\n",T.push(w+"_dtor"),C.push(t[p].Y));xb&&(v+="var ret = retType.fromWireType(rv);\nreturn ret;\n");T.push(v+"}\n");t=Hb(T).apply(null,C);void 0===u[b].X?(t.fa=c-2,u[b]=t):u[b].X[c-2]=t;return[]});return[]})},a:function(a,b,c,d,e,g,k,h,l,m){b=K(b);e=X(d,e);P([],[a],function(n){n=n[0];var r=n.name+"."+b,u={get:function(){Eb("Cannot access "+r+" due to unbound types",[c,k])},enumerable:!0,configurable:!0};l?u.set=function(){Eb("Cannot access "+r+" due to unbound types",[c,k])}:u.set=function(){O(r+
" is a read-only property")};Object.defineProperty(n.S.Z,b,u);P([],l?[c,k]:[c],function(x){var t=x[0],v={get:function(){var C=Ib(this,n,r+" getter");return t.fromWireType(e(g,C))},enumerable:!0};if(l){l=X(h,l);var p=x[1];v.set=function(C){var w=Ib(this,n,r+" setter"),I=[];l(m,w,p.toWireType(I,C));Gb(I)}}Object.defineProperty(n.S.Z,b,v);return[]});return[]})},y:function(a,b){b=K(b);Q(a,{name:b,fromWireType:function(c){var d=Z[c].value;Kb(c);return d},toWireType:function(c,d){return V(d)},argPackAdvance:8,
readValueFromPointer:sb,Y:null})},n:function(a,b,c){c=Ua(c);b=K(b);Q(a,{name:b,fromWireType:function(d){return d},toWireType:function(d,e){if("number"!==typeof e&&"boolean"!==typeof e)throw new TypeError('Cannot convert "'+U(e)+'" to '+this.name);return e},argPackAdvance:8,readValueFromPointer:Lb(b,c),Y:null})},c:function(a,b,c,d,e){function g(m){return m}b=K(b);-1===e&&(e=4294967295);var k=Ua(c);if(0===d){var h=32-8*c;g=function(m){return m<<h>>>h}}var l=-1!=b.indexOf("unsigned");Q(a,{name:b,fromWireType:g,
toWireType:function(m,n){if("number"!==typeof n&&"boolean"!==typeof n)throw new TypeError('Cannot convert "'+U(n)+'" to '+this.name);if(n<d||n>e)throw new TypeError('Passing a number "'+U(n)+'" from JS side to C/C++ side to an argument of type "'+b+'", which is outside the valid range ['+d+", "+e+"]!");return l?n>>>0:n|0},argPackAdvance:8,readValueFromPointer:Mb(b,k,0!==d),Y:null})},b:function(a,b,c){function d(g){g>>=2;var k=G;return new e(za,k[g+1],k[g])}var e=[Int8Array,Uint8Array,Int16Array,Uint16Array,
Int32Array,Uint32Array,Float32Array,Float64Array][b];c=K(c);Q(a,{name:c,fromWireType:d,argPackAdvance:8,readValueFromPointer:d},{sa:!0})},o:function(a,b){b=K(b);var c="std::string"===b;Q(a,{name:b,fromWireType:function(d){var e=G[d>>2];if(c)for(var g=d+4,k=0;k<=e;++k){var h=d+4+k;if(k==e||0==D[h]){g=g?pa(D,g,h-g):"";if(void 0===l)var l=g;else l+=String.fromCharCode(0),l+=g;g=h+1}}else{l=Array(e);for(k=0;k<e;++k)l[k]=String.fromCharCode(D[d+4+k]);l=l.join("")}Y(d);return l},toWireType:function(d,e){e instanceof
ArrayBuffer&&(e=new Uint8Array(e));var g="string"===typeof e;g||e instanceof Uint8Array||e instanceof Uint8ClampedArray||e instanceof Int8Array||O("Cannot pass non-string to std::string");var k=(c&&g?function(){for(var m=0,n=0;n<e.length;++n){var r=e.charCodeAt(n);55296<=r&&57343>=r&&(r=65536+((r&1023)<<10)|e.charCodeAt(++n)&1023);127>=r?++m:m=2047>=r?m+2:65535>=r?m+3:m+4}return m}:function(){return e.length})(),h=Yb(4+k+1);G[h>>2]=k;if(c&&g)qa(e,h+4,k+1);else if(g)for(g=0;g<k;++g){var l=e.charCodeAt(g);
255<l&&(Y(h),O("String has UTF-16 code units that do not fit in 8 bits"));D[h+4+g]=l}else for(g=0;g<k;++g)D[h+4+g]=e[g];null!==d&&d.push(Y,h);return h},argPackAdvance:8,readValueFromPointer:sb,Y:function(d){Y(d)}})},h:function(a,b,c){c=K(c);if(2===b){var d=sa;var e=ua;var g=va;var k=function(){return ta};var h=1}else 4===b&&(d=wa,e=xa,g=ya,k=function(){return G},h=2);Q(a,{name:c,fromWireType:function(l){for(var m=G[l>>2],n=k(),r,u=l+4,x=0;x<=m;++x){var t=l+4+x*b;if(x==m||0==n[t>>h])u=d(u,t-u),void 0===
r?r=u:(r+=String.fromCharCode(0),r+=u),u=t+b}Y(l);return r},toWireType:function(l,m){"string"!==typeof m&&O("Cannot pass non-string to C++ string type "+c);var n=g(m),r=Yb(4+n+b);G[r>>2]=n>>h;e(m,r+4,n+b);null!==l&&l.push(Y,r);return r},argPackAdvance:8,readValueFromPointer:sb,Y:function(l){Y(l)}})},A:function(a,b){b=K(b);Q(a,{ua:!0,name:b,argPackAdvance:0,fromWireType:function(){},toWireType:function(){}})},C:function(a,b,c){a=Nb(a);b=Ob(b,"emval::as");var d=[],e=V(d);F[c>>2]=e;return b.toWireType(d,
a)},k:function(a,b,c,d){a=Rb[a];b=Nb(b);c=Qb(c);a(b,c,null,d)},g:Kb,l:function(a,b){b=Tb(a,b);for(var c=b[0],d=c.name+"_$"+b.slice(1).map(function(m){return m.name}).join("_")+"$",e=["retType"],g=[c],k="",h=0;h<a-1;++h)k+=(0!==h?", ":"")+"arg"+h,e.push("argType"+h),g.push(b[1+h]);d="return function "+Xa("methodCaller_"+d)+"(handle, name, destructors, args) {\n";var l=0;for(h=0;h<a-1;++h)d+="    var arg"+h+" = argType"+h+".readValueFromPointer(args"+(l?"+"+l:"")+");\n",l+=b[h+1].argPackAdvance;d+=
"    var rv = handle[name]("+k+");\n";for(h=0;h<a-1;++h)b[h+1].deleteObject&&(d+="    argType"+h+".deleteObject(arg"+h+");\n");c.ua||(d+="    return retType.toWireType(destructors, rv);\n");e.push(d+"};\n");a=Hb(e).apply(null,g);return Sb(a)},D:function(a,b){a=Nb(a);b=Nb(b);return V(a[b])},p:function(a){4<a&&(Z[a].ka+=1)},r:function(){return V([])},E:function(a){return V(Qb(a))},B:function(a){Gb(Z[a].value);Kb(a)},d:function(a,b){a=Ob(a,"_emval_take_value");a=a.readValueFromPointer(b);return V(a)},
w:function(){A()},u:function(a,b,c){D.copyWithin(a,b,b+c)},v:function(a){a>>>=0;var b=D.length;if(2147483648<a)return!1;for(var c=1;4>=c;c*=2){var d=b*(1+.2/c);d=Math.min(d,a+100663296);d=Math.max(16777216,a,d);0<d%65536&&(d+=65536-d%65536);a:{try{ma.grow(Math.min(2147483648,d)-za.byteLength+65535>>>16);Da();var e=1;break a}catch(g){}e=void 0}if(e)return!0}return!1},x:function(){return 0},s:function(){},m:function(a,b,c,d){for(var e=0,g=0;g<c;g++){for(var k=F[b+8*g>>2],h=F[b+(8*g+4)>>2],l=0;l<h;l++){var m=
D[k+l],n=Ub[a];0===m||10===m?((1===a?ka:B)(pa(n,0)),n.length=0):n.push(m)}e+=h}F[d>>2]=e;return 0},t:function(){},f:function(a){throw Error(a?pa(D,a,void 0):"");}};
(function(){function a(e){f.asm=e.exports;ma=f.asm.G;Da();Ea=f.asm.H;H--;f.monitorRunDependencies&&f.monitorRunDependencies(H);0==H&&(null!==Ka&&(clearInterval(Ka),Ka=null),La&&(e=La,La=null,e()))}function b(e){a(e.instance)}function c(e){return Qa().then(function(g){return WebAssembly.instantiate(g,d)}).then(e,function(g){B("failed to asynchronously prepare wasm: "+g);A(g)})}var d={a:Zb};H++;f.monitorRunDependencies&&f.monitorRunDependencies(H);if(f.instantiateWasm)try{return f.instantiateWasm(d,
a)}catch(e){return B("Module.instantiateWasm callback failed with error: "+e),!1}(function(){return la||"function"!==typeof WebAssembly.instantiateStreaming||Na()||Ma("file://")||"function"!==typeof fetch?c(b):fetch(J,{credentials:"same-origin"}).then(function(e){return WebAssembly.instantiateStreaming(e,d).then(b,function(g){B("wasm streaming compile failed: "+g);B("falling back to ArrayBuffer instantiation");return c(b)})})})();return{}})();
var Xb=f.___wasm_call_ctors=function(){return(Xb=f.___wasm_call_ctors=f.asm.I).apply(null,arguments)},Yb=f._malloc=function(){return(Yb=f._malloc=f.asm.J).apply(null,arguments)},Y=f._free=function(){return(Y=f._free=f.asm.K).apply(null,arguments)};f.___em_js__throwJSError=function(){return(f.___em_js__throwJSError=f.asm.L).apply(null,arguments)};var Db=f.___getTypeName=function(){return(Db=f.___getTypeName=f.asm.M).apply(null,arguments)};
f.___embind_register_native_and_builtin_types=function(){return(f.___embind_register_native_and_builtin_types=f.asm.N).apply(null,arguments)};f.dynCall_jiji=function(){return(f.dynCall_jiji=f.asm.O).apply(null,arguments)};var $b;La=function ac(){$b||bc();$b||(La=ac)};
function bc(){function a(){if(!$b&&($b=!0,f.calledRun=!0,!na)){Ra(Ga);Ra(Ha);if(f.onRuntimeInitialized)f.onRuntimeInitialized();if(f.postRun)for("function"==typeof f.postRun&&(f.postRun=[f.postRun]);f.postRun.length;){var b=f.postRun.shift();Ia.unshift(b)}Ra(Ia)}}if(!(0<H)){if(f.preRun)for("function"==typeof f.preRun&&(f.preRun=[f.preRun]);f.preRun.length;)Ja();Ra(Fa);0<H||(f.setStatus?(f.setStatus("Running..."),setTimeout(function(){setTimeout(function(){f.setStatus("")},1);a()},1)):a())}}
f.run=bc;if(f.preInit)for("function"==typeof f.preInit&&(f.preInit=[f.preInit]);0<f.preInit.length;)f.preInit.pop()();noExitRuntime=!0;bc();