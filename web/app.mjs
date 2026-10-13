import { estimateRt, parseIncidenceCsv } from './rt.mjs';
const get = id => document.getElementById(id);
const ns = 'http://www.w3.org/2000/svg';
let last = null;
const fixed = x => Number.isFinite(x) ? x.toFixed(2) : '—';
function error(message) {
  get('error').textContent = message;
  get('error').hidden = false;
  get('status').textContent = 'Calculation not available';
  get('report').hidden = true;
  last = null;
}
function option(id) {
  const field = get(id);
  if (!field.checkValidity() || field.value.trim() === '') throw new RangeError('Invalid parameter: ' + field.previousElementSibling.textContent);
  return Number(field.value);
}
function svgElement(tag, attrs) {
  const el = document.createElementNS(ns, tag);
  for (const [name, value] of Object.entries(attrs)) el.setAttribute(name, String(value));
  return el;
}
function chart(estimates) {
  const svg = get('chart');
  svg.replaceChildren();
  const left = 56, right = 700, top = 20, bottom = 255;
  const minT = estimates[0].tEnd, maxT = estimates.at(-1).tEnd;
  const maxY = Math.max(1.25, ...estimates.map(r => r.high)) * 1.06;
  const x = t => left + ((t - minT) / Math.max(1, maxT - minT)) * (right - left);
  const y = v => bottom - (v / maxY) * (bottom - top);
  for (let i = 0; i <= 4; i++) {
    const v = i * maxY / 4, yy = y(v);
    svg.append(svgElement('line', {x1:left, x2:right, y1:yy, y2:yy, stroke:'#31435c'}));
    const text = svgElement('text', {x:left-9, y:yy+4,'text-anchor':'end',fill:'#b8c8db','font-size':12});
    text.textContent = fixed(v); svg.append(text);
  }
  const reference = svgElement('line', {x1:left,x2:right,y1:y(1),y2:y(1),stroke:'#f5c777','stroke-width':1.8,'stroke-dasharray':'7 5'});
  svg.append(reference);
  const high = estimates.map((r,i) => (i ? 'L' : 'M') + x(r.tEnd).toFixed(2) + ',' + y(r.high).toFixed(2)).join(' ');
  const low = estimates.slice().reverse().map(r => 'L' + x(r.tEnd).toFixed(2) + ',' + y(r.low).toFixed(2)).join(' ');
  svg.append(svgElement('path', {d:high+' '+low+' Z',fill:'#76bbfa','fill-opacity':'.19',stroke:'none'}));
  svg.append(svgElement('polyline', {points:estimates.map(r => x(r.tEnd)+','+y(r.mean)).join(' '),fill:'none',stroke:'#76bbfa','stroke-width':2.6,'stroke-linejoin':'round'}));
  for (const [xx,txt] of [[left,String(minT)],[(right+left)/2,String(Math.round((minT+maxT)/2))],[right,String(maxT)]]) {
    const t = svgElement('text', {x:xx,y:bottom+23,'text-anchor':'middle',fill:'#b8c8db','font-size':12});
    t.textContent = txt;svg.append(t);
  }
  const axis = svgElement('text',{x:(right+left)/2,y:bottom+40,'text-anchor':'middle',fill:'#b8c8db','font-size':12});
  axis.textContent = 'Day index (zero-based)'; svg.append(axis);
}
function calculate() {
  get('error').hidden = true;
  try {
    const {incidence, dates} = parseIncidenceCsv(get('csv').value);
    const options = {
      siMean:option('siMean'), siSd:option('siSd'), window:option('window'),
      priorMean:option('priorMean'), priorSd:option('priorSd'),
      smoothWindow:option('smoothWindow'), smooth:get('smooth').checked,
    };
    if (options.window > incidence.length) throw new RangeError('Estimation window exceeds the number of days');
    const output = estimateRt(incidence, options);
    if (!output.results.length) throw new RangeError('No estimable windows: preceding infectiousness is zero');
    last = {results:output.results,dates};
    const latest = output.results.at(-1);
    get('latest').textContent = fixed(latest.mean);
    get('interval').textContent = fixed(latest.low)+' – '+fixed(latest.high);
    get('days').textContent = String(output.results.length);
    const tbody = get('rows');tbody.replaceChildren();
    const fragment = document.createDocumentFragment();
    for (const r of output.results) {
      const tr = document.createElement('tr');
      for (const val of [dates[r.tEnd],fixed(r.mean),fixed(r.low),fixed(r.high)]) {
        const td = document.createElement('td'); td.textContent = val; tr.append(td);
      }
      fragment.append(tr);
    }
    tbody.append(fragment);
    chart(output.results);
    get('status').textContent = 'Estimated Rt from '+incidence.length+' daily case observations (latest window ends at '+dates[latest.tEnd]+').';
    get('report').hidden = false;
  } catch (e) {
    error(e instanceof Error ? e.message : 'Unexpected estimation error');
  }
}
get('calculate').addEventListener('click',calculate);
get('csvfile').addEventListener('change',async event=>{
  const file = event.target.files?.[0];
  if (!file) return;
  if (file.size > 1000000) return error('CSV file exceeds the 1 MB limit');
  try {get('csv').value = await file.text(); get('status').textContent = 'File loaded locally. Select Calculate Rt.';get('error').hidden=true;}
  catch {error('Could not read the selected CSV file');}
});
get('example').addEventListener('click',()=>{
  const rows=['date,cases'];
  const date = new Date('2026-01-01T00:00:00Z');
  for (let i=0;i<65;i++) {
    const cases=Math.max(0,Math.round(6*Math.exp(Math.min(i,27)*.095)*Math.exp(-Math.max(i-27,0)*.065)+2*Math.sin(i*1.2)));
    rows.push(date.toISOString().slice(0,10)+','+cases);
    date.setUTCDate(date.getUTCDate()+1);
  }
  get('csv').value=rows.join('\n');
  get('status').textContent='Synthetic example loaded. Select Calculate Rt.';
  get('error').hidden=true;get('report').hidden=true;last=null;
});
get('download').addEventListener('click',()=>{
  if (!last) return;
  const lines=['t_start,t_end,date,mean_r,ci_low_95,ci_high_95,sum_incidence,sum_infectiousness'];
  function cell(s) {return '"'+String(s).replaceAll('"','""')+'"';}
  for (const r of last.results) lines.push([
    r.tStart,r.tEnd,cell(last.dates[r.tEnd]),r.mean,r.low,r.high,r.sumIncidence,r.sumInfectiousness,
  ].join(','));
  const url=URL.createObjectURL(new Blob([lines.join('\r\n')+'\r\n'],{type:'text/csv;charset=utf-8'}));
  const a=document.createElement('a');a.href=url;a.download='rt_estimates.csv';a.click();
  setTimeout(()=>URL.revokeObjectURL(url),1000);
});
get('example').click();
