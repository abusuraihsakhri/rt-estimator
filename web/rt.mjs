// Cori et al. (2013) renewal-equation estimator, implemented for offline browser use.
// Gamma CDF, inverse CDF and serial interval follow the Python implementation.
const lanczos = [
  0.99999999999980993, 676.5203681218851, -1259.1392167224028,
  771.32342877765313, -176.61502916214059, 12.507343278686905,
  -0.13857109526572012, 9.984369578019572e-6, 1.5056327351493116e-7,
];

function logGamma(z) {
  if (z < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1 - z);
  let x = lanczos[0];
  z -= 1;
  for (let i = 1; i < lanczos.length; i++) x += lanczos[i] / (z + i);
  const t = z + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

function gammaP(a, x) {
  if (!(a > 0) || !Number.isFinite(x) || x < 0) throw new RangeError('Invalid gamma CDF arguments');
  if (x === 0) return 0;
  const factor = Math.exp(a * Math.log(x) - x - logGamma(a));
  if (x < a + 1) {
    let sum = 1 / a;
    let term = sum;
    for (let n = 1; n <= 10000; n++) {
      term *= x / (a + n);
      sum += term;
      if (Math.abs(term) <= Math.abs(sum) * 1e-14) break;
      if (n === 10000) throw new RangeError('Gamma CDF did not converge');
    }
    return Math.min(1, Math.max(0, factor * sum));
  }
  const tiny = 1e-300;
  let b = x + 1 - a;
  let c = 1 / tiny;
  let d = 1 / (Math.abs(b) < tiny ? tiny : b);
  let h = d;
  for (let i = 1; i <= 10000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < tiny) d = tiny;
    c = b + an / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const delta = d * c;
    h *= delta;
    if (Math.abs(delta - 1) < 1e-14) break;
    if (i === 10000) throw new RangeError('Gamma CDF did not converge');
  }
  return Math.min(1, Math.max(0, 1 - factor * h));
}

export function gammaCdf(x, shape, scale = 1) {
  if (!(shape > 0 && scale > 0) || !Number.isFinite(shape) || !Number.isFinite(scale)) {
    throw new RangeError('Gamma shape and scale must be finite and positive');
  }
  return gammaP(shape, Math.max(0, x) / scale);
}

export function gammaQuantile(prob, shape, scale = 1) {
  if (!(prob > 0 && prob < 1) || !Number.isFinite(prob)) throw new RangeError('Probability must be between 0 and 1');
  if (!(shape > 0 && scale > 0) || !Number.isFinite(shape) || !Number.isFinite(scale)) throw new RangeError('Invalid posterior parameters');
  let low = 0;
  let high = Math.max(1, shape);
  while (gammaP(shape, high) < prob) {
    high *= 2;
    if (!Number.isFinite(high)) throw new RangeError('Posterior quantile overflow');
  }
  for (let i = 0; i < 75; i++) {
    const mid = (low + high) / 2;
    if (gammaP(shape, mid) < prob) low = mid; else high = mid;
  }
  return (low + high) * 0.5 * scale;
}

export function discretizeSerialInterval(mean = 4.7, sd = 2.9, maxDays = null) {
  if (!Number.isFinite(mean) || mean <= 1 || !Number.isFinite(sd) || sd <= 0) {
    throw new RangeError('Serial interval mean must be > 1 and SD must be positive and finite');
  }
  const days = maxDays === null ? Math.ceil(mean + 10 * sd) : maxDays;
  if (!Number.isInteger(days) || days < 1 || days > 365) throw new RangeError('Serial interval support must be 1–365 days');
  const shape = ((mean - 1) / sd) ** 2;
  const scale = sd * sd / (mean - 1);
  const F = (k, a) => k < 0 ? 0 : gammaCdf(k, a, scale);
  const w = [0];
  for (let k = 1; k <= days; k++) {
    let p = k * F(k, shape) + (k - 2) * F(k - 2, shape) - 2 * (k - 1) * F(k - 1, shape);
    p += shape * scale * (2 * F(k - 1, shape + 1) - F(k - 2, shape + 1) - F(k, shape + 1));
    w.push(Math.max(0, p));
  }
  const total = w.reduce((a, b) => a + b, 0);
  if (!Number.isFinite(total) || !(total > 0)) throw new RangeError('Invalid serial interval distribution');
  return w.map(v => v / total);
}

export function trailingAverage(x, window = 7) {
  if (!Number.isInteger(window) || window < 1) throw new RangeError('Smoothing window must be positive');
  const out = [];
  let sum = 0;
  for (let t = 0; t < x.length; t++) {
    sum += x[t];
    if (t >= window) sum -= x[t - window];
    out.push(sum / Math.min(window, t + 1));
  }
  return out;
}

export function estimateRt(incidence, {window = 7, siMean = 4.7, siSd = 2.9, priorMean = 5, priorSd = 5, smooth = false, smoothWindow = 7} = {}) {
  if (!Array.isArray(incidence) || incidence.length < 2 || incidence.length > 5000 ||
      incidence.some(x => !Number.isFinite(x) || x < 0)) throw new RangeError('Provide 2–5000 finite non-negative incidence values');
  if (!Number.isInteger(window) || window < 1 || window > incidence.length) throw new RangeError('Invalid estimation window');
  if (!Number.isFinite(priorMean) || priorMean <= 0 || !Number.isFinite(priorSd) || priorSd <= 0) {
    throw new RangeError('Prior mean and SD must be finite and positive');
  }
  const si = discretizeSerialInterval(siMean, siSd);
  const data = smooth ? trailingAverage(incidence, smoothWindow) : incidence.slice();
  const lam = new Array(data.length).fill(0);
  for (let t = 1; t < data.length; t++) {
    for (let s = 1; s < si.length && s <= t; s++) lam[t] += si[s] * data[t - s];
  }
  const priorShape = (priorMean / priorSd) ** 2;
  const priorRate = priorMean / (priorSd * priorSd);
  const results = [];
  let sumI = 0, sumLam = 0;
  for (let t = 0; t < data.length; t++) {
    sumI += data[t]; sumLam += lam[t];
    if (t >= window) { sumI -= data[t - window]; sumLam -= lam[t - window]; }
    if (t < window - 1 || sumLam <= 0) continue;
    const shape = priorShape + sumI;
    const scale = 1 / (priorRate + sumLam);
    results.push({
      tStart: t - window + 1, tEnd: t, mean: shape * scale,
      low: gammaQuantile(0.025, shape, scale), high: gammaQuantile(0.975, shape, scale),
      sumIncidence: sumI, sumInfectiousness: sumLam,
    });
  }
  return {results, incidence: data, infectiousness: lam, serialInterval: si};
}

// Parse RFC 4180 fields including quoted newlines, commas, and escaped quotes.
export function parseIncidenceCsv(text) {
  if (typeof text !== 'string' || text.length > 1000000) throw new RangeError('CSV must be text under 1 MB');
  text = text.replace(/^\uFEFF/, '');
  const rows = []; let row = [], field = '', quoted = false, closedQuote = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (ch === '"') { quoted = false; closedQuote = true; }
      else field += ch;
    } else if (ch === '"' && field === '' && !closedQuote) quoted = true;
    else if (ch === ',' || ch === '\r' || ch === '\n') {
      row.push(field); field = ''; closedQuote = false;
      if (ch !== ',') {
        if (row.some(v => v.trim() !== '')) rows.push(row);
        row = [];
        if (ch === '\r' && text[i + 1] === '\n') i++;
      }
    } else {
      if (closedQuote && ch !== ' ' && ch !== '\t') throw new RangeError('Unexpected text after closing CSV quote');
      field += ch;
    }
  }
  if (quoted) throw new RangeError('Unclosed CSV quote');
  row.push(field);
  if (row.some(v => v.trim() !== '')) rows.push(row);
  if (rows.length < 3) throw new RangeError('CSV must include a header and at least two daily values');
  const names = rows[0].map(x => x.trim().toLowerCase());
  const caseIndex = ['cases', 'incidence', 'count'].map(n => names.indexOf(n)).find(i => i >= 0);
  if (caseIndex === undefined) throw new RangeError('CSV header must include cases, incidence or count');
  const dateIndex = names.indexOf('date');
  if (rows.length - 1 > 5000) throw new RangeError('Maximum 5000 days per analysis');
  const incidence = [], dates = [];
  for (let i = 1; i < rows.length; i++) {
    const cells = rows[i];
    if (cells.length !== names.length) throw new RangeError('Row ' + (i + 1) + ' has an incorrect number of columns');
    const raw = cells[caseIndex].trim();
    const n = Number(raw);
    if (!raw || !Number.isSafeInteger(n) || n < 0 || n > 10000000) {
      throw new RangeError('Row ' + (i + 1) + ': cases must be non-negative integers up to 10 million');
    }
    incidence.push(n); dates.push(dateIndex >= 0 ? cells[dateIndex].trim() : String(i - 1));
  }
  return {incidence, dates};
}
