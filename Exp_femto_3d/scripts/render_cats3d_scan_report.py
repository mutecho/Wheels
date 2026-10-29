#!/usr/bin/env python3
"""Build an offline browser index over ROOT-rendered, saved CATS3D scan points."""

import argparse
import json
import struct
from pathlib import Path


def png_payload(path):
    """Reuse ROOT's lossless PNG IDAT stream as a PDF Flate image (PNG predictor).

    ROOT emits non-interlaced 8-bit RGB/gray images. No pixel resampling or extra
    image/PDF package is needed, and the original vector canvases remain in ROOT.
    """
    raw = path.read_bytes()
    if raw[:8] != b'\x89PNG\r\n\x1a\n':
        raise ValueError('Invalid PNG: ' + str(path))
    cursor, payload, dimensions = 8, [], None
    while cursor + 12 <= len(raw):
        size = struct.unpack_from('>I', raw, cursor)[0]
        kind = raw[cursor + 4:cursor + 8]
        value = raw[cursor + 8:cursor + 8 + size]
        if kind == b'IHDR':
            width, height, depth, color, compression, filtering, interlace = struct.unpack('>IIBBBBB', value)
            if depth != 8 or color not in (0, 2) or compression or filtering or interlace:
                raise ValueError('Unsupported ROOT PNG encoding: ' + str(path))
            dimensions = width, height, 1 if color == 0 else 3
        elif kind == b'IDAT':
            payload.append(value)
        elif kind == b'IEND':
            break
        cursor += size + 12
    if dimensions is None or not payload:
        raise ValueError('Incomplete PNG: ' + str(path))
    return (*dimensions, b''.join(payload))


def write_report_pdf(directory, report, points):
    """Place the ROOT images on complete landscape pages; keep all failed points."""
    objects = [None, None, b'<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>']
    pages = []

    def add(value):
        objects.append(value)
        return len(objects)

    def stream(value, attributes=b''):
        return b'<< /Length %d ' % len(value) + attributes + b' >>\nstream\n' + value + b'\nendstream'

    def page(content, resources):
        content_id = add(stream(content))
        pages.append(add((f'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 842 595] '
                          f'/Resources << {resources} >> /Contents {content_id} 0 R >>').encode()))

    title = ('fixed final lambda' if report['kind'] == 'constant' else 'fixed input-table scale')
    lines = [
        'CATS3D Gaussian scan report',
        f'{title} = {report["coordinate"]}; alpha = 2',
        f'{report["points"]} saved slice points; {report["valid"]} valid minima; '
        f'{report["errors_valid"]} valid conditional errors.',
        'Orange solid: valid saved fit. Gray dashed: unconverged diagnostic trial.',
        'Saved parameters, statuses and all seed attempts are retained. No minimization is rerun.',
        'Data/model use the same transverse normal bins: |q_other| < 0.04 GeV/c.',
        'Model curves cover computed q bins only. No extra rebinning or uncertainty band.',
        'Errors condition on fixed external lambda and the numerical kernel.',
        'No automatic production lambda selection or confidence-interval interpretation.',
        'Group parameter pages precede the three-projection page for every slice.',
        'PDF: lossless ROOT image snapshots. Original vector canvases are in the companion ROOT file.',
    ]
    text = []
    for i, line in enumerate(lines):
        escaped = line.replace('\\', '\\\\').replace('(', r'\(').replace(')', r'\)')
        text.append(f'BT /F1 {25 if i == 0 else 13} Tf 45 {535 - i * 36} Td ({escaped}) Tj ET')
    page('\n'.join(text).encode('ascii'), '/Font << /F1 3 0 R >>')
    images = report['group_images'] + [p['image'] for p in points]
    for name in images:
        width, height, colors, payload = png_payload(directory / name)
        space = '/DeviceGray' if colors == 1 else '/DeviceRGB'
        attributes = (f'/Type /XObject /Subtype /Image /Width {width} /Height {height} '
                      f'/ColorSpace {space} /BitsPerComponent 8 /Filter /FlateDecode '
                      f'/DecodeParms << /Predictor 15 /Colors {colors} /BitsPerComponent 8 /Columns {width} >>').encode()
        image_id = add(stream(payload, attributes))
        scale = min(810 / width, 563 / height)
        w, h = width * scale, height * scale
        x, y = (842 - w) / 2, (595 - h) / 2
        page(f'q {w:.6f} 0 0 {h:.6f} {x:.6f} {y:.6f} cm /Im0 Do Q'.encode(),
             f'/XObject << /Im0 {image_id} 0 R >>')
    if len(pages) != report['pages']:
        raise ValueError('Report page coverage mismatch: ' + report['stem'])
    objects[0] = b'<< /Type /Catalog /Pages 2 0 R >>'
    objects[1] = (f'<< /Type /Pages /Count {len(pages)} /Kids [ '
                  + ' '.join(f'{p} 0 R' for p in pages) + ' ] >>').encode()
    target = directory / (report['stem'] + '.pdf')
    temporary = target.with_suffix('.pdf.tmp')
    with temporary.open('wb') as output:
        output.write(b'%PDF-1.4\n%\xe2\xe3\xcf\xd3\n')
        offsets = [0]
        for i, value in enumerate(objects, 1):
            offsets.append(output.tell())
            output.write(f'{i} 0 obj\n'.encode() + value + b'\nendobj\n')
        startxref = output.tell()
        output.write(f'xref\n0 {len(offsets)}\n0000000000 65535 f \n'.encode())
        for offset in offsets[1:]:
            output.write(f'{offset:010d} 00000 n \n'.encode())
        output.write((f'trailer\n<< /Size {len(offsets)} /Root 1 0 R >>\n'
                      f'startxref\n{startxref}\n%%EOF\n').encode())
    temporary.replace(target)


TEMPLATE = r'''<!doctype html>
<html lang="zh-CN">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>CATS3D · λ 扫描拟合报表</title>
<style>
:root{--bg:#f5f7f2;--paper:#fffef9;--soft:#eef5ef;--ink:#1e2621;--muted:#647067;--line:#d9e1d8;--accent:#28786a;--warn:#986037;--bad:#a74038}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.7 ui-sans-serif,system-ui,-apple-system,"PingFang SC",sans-serif}
main{max-width:1100px;margin:auto;padding:32px 28px 64px}a{color:var(--accent);text-underline-offset:3px}h1{font-size:28px;margin:4px 0 12px}h2{font-size:21px;margin:30px 0 12px}h3{font-size:17px;margin:0}p{max-width:80ch;margin:8px 0 16px}.eyebrow{font-size:12px;color:var(--accent);font-weight:650}code{font:12px/1.5 ui-monospace,monospace;background:var(--soft);padding:2px 5px;border-radius:4px;overflow-wrap:anywhere}.muted{color:var(--muted)}
.downloads{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:22px 0}.download{border:1px solid var(--line);border-radius:8px;background:var(--paper);padding:16px}.download p{font-size:13px;color:var(--muted);margin:8px 0}.download a{margin-right:14px}.filters{display:flex;gap:12px;flex-wrap:wrap;padding:16px;border:1px solid var(--line);border-radius:8px;background:var(--paper)}label{font-size:12px;color:var(--muted);display:flex;flex-direction:column;gap:4px}select{font:inherit;font-size:14px;padding:8px;border:1px solid var(--line);border-radius:5px;background:var(--paper);color:var(--ink);max-width:400px}
.table-wrap{overflow:auto;max-height:420px;border:1px solid var(--line);border-radius:7px;background:var(--paper);margin-top:14px}table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}th,td{text-align:left;border-bottom:1px solid var(--line);padding:9px 12px;white-space:nowrap}th{position:sticky;top:0;background:var(--soft);font-weight:600}tr.pick{cursor:pointer}tr.pick:hover,tr.selected{background:var(--soft)}.bad{color:var(--bad)}.warn{color:var(--warn)}.ok{color:var(--accent)}
.viewer{margin:18px 0;border:1px solid var(--line);border-radius:8px;overflow:hidden;background:white}.viewer img{width:100%;display:block}.viewer figcaption{padding:12px 18px;border-top:1px solid var(--line);font-size:13px}.detail{display:grid;grid-template-columns:1fr 1fr;gap:18px;align-items:start}.detail .table-wrap{margin:0;max-height:none}.status-note{border-left:3px solid var(--accent);padding:8px 16px;background:var(--soft);margin:14px 0}.status-note.bad{border-color:var(--bad);background:#faf0ed}footer{margin-top:32px;border-top:1px solid var(--line);padding-top:18px;font-size:12px;color:var(--muted)}footer p{max-width:none}.source{font-size:12px;overflow-wrap:anywhere}.count{font-size:13px;margin-top:12px}#empty{padding:24px;color:var(--muted)}button{border:1px solid var(--line);background:var(--paper);border-radius:4px;padding:7px 11px;color:var(--accent);cursor:pointer;font:inherit}
@media(max-width:760px){main{padding:22px 14px}.downloads{grid-template-columns:repeat(2,1fr)}h1{font-size:23px}.detail{grid-template-columns:1fr}select{max-width:calc(100vw - 64px)}.filters{gap:10px}}
@media print{.filters,.table-wrap,#empty{display:none}.viewer{break-inside:avoid}.downloads{grid-template-columns:repeat(4,1fr)}body{background:white}}
</style>
<main>
<header><div class="eyebrow">EXP_FEMTO_3D / CATS3D / FIXED EXTERNAL λ</div>
<h1>λ 扫描：逐点拟合报表</h1>
<p>高斯源 α = 2。按保存参数重建三维模型与 out / side / long 投影，保留扫描原有的收敛状态和条件误差。每个 PDF 含组内参数随 φ 的分布，以及全部切片的拟合投影。</p>
<p class="muted" id="overview"></p></header>
<nav class="downloads" id="downloads" aria-label="各扫描坐标的报表"></nav>
<h2>查看切片</h2>
<div class="filters">
<label>扫描坐标<select id="coordinate"></select></label>
<label>中心度 / mT 组<select id="group"></select></label>
<label>结果状态<select id="status"><option value="all">全部状态</option><option value="valid">有效解</option><option value="errors">有效解且误差有效</option><option value="noerrors">有效解但误差不可用</option><option value="failed">未获得有效解</option></select></label>
</div>
<div class="count muted" id="count"></div>
<div class="table-wrap"><table><thead><tr><th>中心度 / mT</th><th>φ 切片</th><th>解 / 误差</th><th>目标函数 / ndf</th><th>Δ目标函数</th><th>EDM</th><th>尝试次数</th></tr></thead><tbody id="rows"></tbody></table></div>
<div id="empty" hidden>当前筛选条件没有切片。</div>
<section id="selected">
<div class="status-note" id="status-note"></div>
<figure class="viewer"><img id="image" alt="三个动量方向的 CF 数据与保存模型投影，以及拟合参数和状态"><figcaption id="caption"></figcaption></figure>
<div class="detail"><div><h3>保存参数</h3><p class="muted">半径矩阵单位为 fm²，q² baseline 系数单位为 GeV⁻²。固定输入不计入自由参数。</p><div class="table-wrap"><table><thead><tr><th>参数</th><th>值</th><th>条件误差 / 状态</th></tr></thead><tbody id="parameters"></tbody></table></div></div>
<div><h3>诊断与文件</h3><p id="diagnostics"></p><p id="links"></p><p class="muted">橙色实线表示有效的保存拟合解；灰色虚线表示未收敛试探结果，仅供诊断。误差缺失保留为空，不补零；投影未绘制参数误差带。</p><p class="muted">数据与模型平均同一组横向 bins（|q| &lt; 0.04 GeV/c），模型曲线只覆盖已计算的原始 q bins。这里没有额外重分箱。</p></div></div>
</section>
<footer><p>扫描结果不自动选择正式 λ，也不把 Δ目标函数标为置信区间。HESSE 误差条件于固定的外部 λ 和数值核。组目标函数要求全部非积分 φ 成员有效，φ 积分切片单独保存。</p>
<p><a href="scan_points.csv" download>全部参数 CSV</a> · <a href="manifest.json" download>报表索引与状态 JSON</a></p>
<p class="source" id="source"></p><p id="reconstruction"></p></footer>
</main>
<script id="report-data" type="application/json">__REPORT_DATA__</script>
<script>
const data=JSON.parse(document.getElementById('report-data').textContent);
const byId=id=>document.getElementById(id),fmt=x=>x===null?'—':String(Number(Number(x).toPrecision(5)));
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const key=p=>p.kind+'|'+p.coordinate;
const label=p=>p.kind==='constant'?'λ = '+p.coordinate:'表倍率 s = '+p.coordinate;
const groupLabel=s=>s.replace(/^cent_/,'cent ').replace(/__mtbin_\d+__mt_/,' / mT ').replace(/__/g,' / ');
const reports=new Map(data.reports.map(r=>[key(r),r]));
let selected=null;
const valid=data.points.filter(p=>p.valid).length,errors=data.points.filter(p=>p.valid&&p.errors_valid).length;
byId('overview').textContent=`${data.reports.length} 个扫描坐标 · ${data.points.length} 个切片点 · ${valid} 个有效解 · ${errors} 个条件误差有效。全部失败记录保留。`;
byId('downloads').innerHTML=data.reports.map(r=>`<div class="download"><h3>${escape(label(r))}</h3><p>${r.valid} / ${r.points} 有效解<br>${r.errors_valid} 个误差有效 · ${r.pages} 页</p><a href="${r.stem}.pdf">PDF</a><a href="${r.stem}.root" download>ROOT</a></div>`).join('');
for(const r of data.reports)byId('coordinate').add(new Option(label(r),key(r)));
const groupIds=[...new Set(data.points.map(p=>p.group))];
byId('group').add(new Option('全部组','all'));for(const g of groupIds)byId('group').add(new Option(groupLabel(g),g));
const initial=new URLSearchParams(location.hash.slice(1));
if(reports.has(initial.get('coordinate')))byId('coordinate').value=initial.get('coordinate');
if(groupIds.includes(initial.get('group')))byId('group').value=initial.get('group');
function show(p){
 selected=p;
 const r=reports.get(key(p));
 byId('image').src=p.image;
 byId('caption').textContent=label(p)+' · '+groupLabel(p.group)+' · '+(p.phi_integrated?'φ 积分':'φ = '+fmt(p.phi)+' rad');
 const note=byId('status-note');note.className='status-note'+(p.valid?'':' bad');
 note.textContent=p.valid?(p.errors_valid?'最小化与条件误差均有效。':'最小化有效；条件误差不可用。'):'未获得有效解：'+p.status+'。图中的试探曲线不能视为有效拟合结果。';
 byId('parameters').innerHTML=p.parameters.filter(v=>!v.fixed||v.name==='alpha'||v.name.includes('Alpha')||(v.name.includes('Carrier')&&p.kind==='constant')).map(v=>`<tr><td>${escape(v.name)}</td><td>${fmt(v.value)}</td><td class="${v.at_bound?'warn':''}">${v.fixed?'固定':p.valid&&p.errors_valid?fmt(v.error):'不可用'}${v.at_bound?' · 触及边界':''}</td></tr>`).join('');
 byId('diagnostics').textContent=`MIGRAD ${p.migrad} · HESSE ${p.hesse} · covariance ${p.covariance} · EDM ${fmt(p.edm)}。保存目标函数 ${fmt(p.objective)}，ndf ${p.ndf}。${p.attempts} 次种子尝试。${p.reason?' '+p.reason:''}${p.reconstruction_error?' 重建：'+p.reconstruction_error:''}`;
 byId('links').innerHTML=`<a href="${r.stem}.pdf">本 λ 的 PDF</a> · <a href="${r.stem}.root" download>ROOT 报表</a> · <a href="${p.image}">原始图片</a>`;
 document.querySelectorAll('tr.pick').forEach(row=>row.classList.toggle('selected',row.dataset.id===p.id));
}
function update(){
 const k=byId('coordinate').value,g=byId('group').value,status=byId('status').value;
 const points=data.points.filter(p=>key(p)===k&&(g==='all'||p.group===g)&&
  (status==='all'||status==='valid'&&p.valid||status==='errors'&&p.valid&&p.errors_valid||status==='noerrors'&&p.valid&&!p.errors_valid||status==='failed'&&!p.valid));
 byId('count').textContent=`显示 ${points.length} 个切片点，其中 ${points.filter(p=>p.valid).length} 个有效。`;
 byId('rows').innerHTML=points.map((p,i)=>`<tr class="pick" data-index="${i}" data-id="${escape(p.id)}"><td>${escape(groupLabel(p.group))}</td><td>${p.phi_integrated?'积分':fmt(p.phi)+' rad'}</td><td class="${p.valid?p.errors_valid?'ok':'warn':'bad'}">${p.valid?'有效':'失败'} / ${p.valid&&p.errors_valid?'有效':'不可用'}</td><td>${p.valid&&p.ndf>0?fmt(p.objective/p.ndf):'—'}</td><td>${fmt(p.delta_objective)}</td><td>${fmt(p.edm)}</td><td>${p.attempts}</td></tr>`).join('');
 byId('rows').querySelectorAll('tr').forEach(row=>row.addEventListener('click',()=>show(points[Number(row.dataset.index)])));
 byId('empty').hidden=points.length>0;byId('selected').hidden=points.length===0;
 if(points.length)show(points.find(p=>selected&&p.id===selected.id)||points.find(p=>p.valid&&p.errors_valid)||points[0]);
 history.replaceState(null,'','#'+new URLSearchParams({coordinate:k,group:g}));
}
for(const id of ['coordinate','group','status'])byId(id).addEventListener('change',update);
byId('source').textContent='扫描输入：'+data.scan_root+'\n原始 CF：'+data.input_cf_root;
byId('reconstruction').textContent=`本次报表重建用时 ${Math.round(data.elapsed_seconds)} 秒，进行了 ${data.source_folds} 次源卷积；没有运行 Minuit 或改写扫描结果。模型使用扫描文件记录的原生 kernel。`;
update();
</script>
</html>
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report_directory', type=Path)
    args = parser.parse_args()
    manifest = args.report_directory / 'manifest.json'
    data = json.loads(manifest.read_text())
    if not data.get('complete'):
        raise SystemExit('Report reconstruction is not complete.')
    for report in data['reports']:
        if not (args.report_directory / (report['stem'] + '.root')).is_file():
            raise SystemExit('Missing report ROOT: ' + report['stem'])
        points = [p for p in data['points'] if p['kind'] == report['kind'] and p['coordinate'] == report['coordinate']]
        write_report_pdf(args.report_directory, report, points)
    for point in data['points']:
        if not (args.report_directory / point['image']).is_file():
            raise SystemExit('Missing point image: ' + point['image'])
    encoded = json.dumps(data, ensure_ascii=False, allow_nan=False).replace('<', r'\u003c')
    target = args.report_directory / 'index.html'
    target.write_text(TEMPLATE.replace('__REPORT_DATA__', encoded))
    print(target)


if __name__ == '__main__':
    main()
