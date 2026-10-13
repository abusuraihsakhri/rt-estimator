import test from 'node:test';
import assert from 'node:assert/strict';
import { gammaCdf, gammaQuantile, discretizeSerialInterval, estimateRt, parseIncidenceCsv, trailingAverage } from './rt.mjs';

test('gamma CDF and inverse agree with analytic exponential and known values', () => {
  assert.ok(Math.abs(gammaCdf(1, 1) - (1 - Math.exp(-1))) < 1e-11);
  assert.ok(Math.abs(gammaQuantile(0.5, 1) - Math.log(2)) < 1e-10);
  assert.ok(Math.abs(gammaQuantile(0.975, 2, 3) - 16.7149301728) < 0.00001);
});

test('serial interval discretization is normalized and agrees with Python implementation', () => {
  const si = discretizeSerialInterval(4.7, 2.9);
  assert.equal(si[0], 0);
  assert.ok(si.every(v => v >= 0));
  assert.ok(Math.abs(si.reduce((a,b)=>a+b,0) - 1) < 1e-12);
  const mean = si.reduce((a,v,k) => a + k*v, 0);
  assert.ok(Math.abs(mean - 4.7) < 0.5);
});

test('constant incidence yields R near 1 after transient', () => {
  const {results} = estimateRt(new Array(50).fill(100));
  assert.ok(results.length > 0);
  assert.ok(Math.abs(results.at(-1).mean - 1) < 0.02);
  assert.ok(results.at(-1).low < 1 && results.at(-1).high > 1);
});

test('zero incidence and invalid inputs fail safely', () => {
  assert.equal(estimateRt(new Array(20).fill(0)).results.length, 0);
  assert.throws(() => estimateRt([1, -1, 3]));
  assert.throws(() => estimateRt([1, NaN, 3]));
  assert.throws(() => estimateRt([1, 2], {window: 3}));
  assert.throws(() => discretizeSerialInterval(Infinity, 2));
});

test('causal smoothing and CSV handling', () => {
  assert.deepEqual(trailingAverage([0, 10, 20], 2), [0, 5, 15]);
  const x = parseIncidenceCsv('date,cases\r\n"2026-01-01",10\r\n"2026-01-02",12\r\n');
  assert.deepEqual(x.incidence, [10, 12]);
  assert.deepEqual(x.dates, ['2026-01-01', '2026-01-02']);
  assert.deepEqual(parseIncidenceCsv('count,comment\n1,"a,b"\n2,"c,d"').incidence, [1,2]);
  assert.throws(() => parseIncidenceCsv('cases\n1\nNaN'));
  assert.throws(() => parseIncidenceCsv('cases\n1\n'));
  assert.throws(() => parseIncidenceCsv('cases\n1\n-2'));
});

test('browser posterior matches independent SciPy reference values', () => {
  const si = discretizeSerialInterval(4.7, 2.9);
  assert.ok(Math.abs(si[1] - 0.05650093) < 2e-7);
  assert.ok(Math.abs(si[2] - 0.17807472) < 2e-7);
  const x = [10,12,11,13,14,20,18,21,23,22,27,30,28,25,29,31,33,35,34,35];
  const final = estimateRt(x).results.at(-1);
  assert.ok(Math.abs(final.mean - 1.1991399278771875) < 2e-6);
  assert.ok(Math.abs(final.low - 1.0469074916220034) < 2e-6);
  assert.ok(Math.abs(final.high - 1.3615553428548552) < 2e-6);
});
