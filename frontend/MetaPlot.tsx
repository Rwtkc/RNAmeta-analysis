import * as Select from '@radix-ui/react-select'
import { Check, ChevronDown, Download, LineChart, SlidersHorizontal } from 'lucide-react'
import { jsPDF } from 'jspdf'
import { useEffect, useMemo, useRef, useState } from 'react'
import { svg2pdf } from 'svg2pdf.js'
import type { ParsedBed } from './rnaenrich-core'

const appBase = import.meta.env.BASE_URL
const apiBase = `${appBase}RNAenrich/api`

type SpeciesOption = { value: string; label: string }
type MetaPlotControls = {
  pltTxType: 'tx' | 'mrna' | 'ncrna'
  headOrtail: 'TRUE' | 'FALSE'
  txPrimaryOnly: 'TRUE' | 'FALSE'
  txlncrnaOverlapmrna: 'TRUE' | 'FALSE'
  mapFilterTranscript: 'TRUE' | 'FALSE'
  enableCI: 'TRUE' | 'FALSE'
  stSampleModle: 'Equidistance' | 'random'
  stSampleNum: number
  stAmblguity: number
  txfiveutrMinLength: number
  txcdsMinLength: number
  txthreeutrMinLength: number
  txlongNcrnaMinLength: number
  txpromoterLength: number
  txtailLength: number
  adjust: number
  overlapIndex: number
  siteLengthIndex: number
}
type ChartPoint = { x: number; density: number; confidenceDown?: number | null; confidenceUp?: number | null }
type ChartSeries = { name: string; originalName: string; color: string; values: ChartPoint[] }
type ChartSegment = { start: number; end: number; mid: number; label: string; component: string; alpha: number; height: number }
type ChartPayload = {
  title: string
  yLabel: string
  showCI: boolean
  xDomain?: [number, number]
  yDomain?: [number, number]
  components: { trackBase?: number; labelY?: number; segments: ChartSegment[]; separators: { x: number }[] }
  series: ChartSeries[]
}
type MetaPlotResult = {
  chart: ChartPayload
  summary: { species: string; sample_count: number; interval_count: number; sampled_point_count: number; transcript_hits: number }
  data: Record<string, string | number | null>[]
}

const defaults: MetaPlotControls = {
  pltTxType: 'tx', headOrtail: 'TRUE', txPrimaryOnly: 'TRUE', txlncrnaOverlapmrna: 'TRUE',
  mapFilterTranscript: 'TRUE', enableCI: 'FALSE', stSampleModle: 'Equidistance', stSampleNum: 10,
  stAmblguity: 5, txfiveutrMinLength: 100, txcdsMinLength: 100, txthreeutrMinLength: 100,
  txlongNcrnaMinLength: 100, txpromoterLength: 1000, txtailLength: 1000, adjust: 1,
  overlapIndex: 1, siteLengthIndex: 1,
}

const chartSize = { width: 1180, height: 520 }
const chartMargin = { top: 58, right: 28, bottom: 112, left: 76 }

type Point = { x: number; y: number }
type ChartHover = { series: ChartSeries; point: ChartPoint; left: number; top: number }

function basisPath(points: Point[]) {
  if (!points.length) return ''
  if (points.length === 1) return `M${points[0].x},${points[0].y}`
  if (points.length === 2) return `M${points[0].x},${points[0].y}L${points[1].x},${points[1].y}`
  const format = (value: number) => value.toFixed(2)
  let previous = points[0]
  let current = points[1]
  let path = `M${format(previous.x)},${format(previous.y)}L${format((5 * previous.x + current.x) / 6)},${format((5 * previous.y + current.y) / 6)}`
  for (let index = 2; index < points.length; index += 1) {
    const next = points[index]
    path += `C${format((2 * previous.x + current.x) / 3)},${format((2 * previous.y + current.y) / 3)},${format((previous.x + 2 * current.x) / 3)},${format((previous.y + 2 * current.y) / 3)},${format((previous.x + 4 * current.x + next.x) / 6)},${format((previous.y + 4 * current.y + next.y) / 6)}`
    previous = current
    current = next
  }
  path += `C${format((2 * previous.x + current.x) / 3)},${format((2 * previous.y + current.y) / 3)},${format((previous.x + 2 * current.x) / 3)},${format((previous.y + 2 * current.y) / 3)},${format((previous.x + 5 * current.x) / 6)},${format((previous.y + 5 * current.y) / 6)}L${format(current.x)},${format(current.y)}`
  return path
}

function chartPoints(values: ChartPoint[], xScale: (value: number) => number, yScale: (value: number) => number, key: 'density' | 'confidenceDown' | 'confidenceUp') {
  return values.flatMap((point) => {
    const value = point[key]
    return value === null || value === undefined || !Number.isFinite(value) ? [] : [{ x: xScale(point.x), y: yScale(value) }]
  })
}

function linePath(values: ChartPoint[], xScale: (value: number) => number, yScale: (value: number) => number, key: 'density' | 'confidenceDown' | 'confidenceUp') {
  return basisPath(chartPoints(values, xScale, yScale, key))
}

function areaPath(values: ChartPoint[], xScale: (value: number) => number, yScale: (value: number) => number) {
  const points = chartPoints(values, xScale, yScale, 'density')
  if (!points.length) return ''
  return `${basisPath(points)}L${points.at(-1)?.x},${yScale(0)}L${points[0].x},${yScale(0)}Z`
}

function confidenceAreaPath(values: ChartPoint[], xScale: (value: number) => number, yScale: (value: number) => number) {
  const upper = chartPoints(values, xScale, yScale, 'confidenceUp')
  const lower = chartPoints([...values].reverse(), xScale, yScale, 'confidenceDown')
  if (!upper.length || !lower.length) return ''
  return `${basisPath(upper)}${basisPath(lower).replace(/^M/, 'L')}Z`
}

function tickValues(start: number, stop: number, count = 5) {
  const rawStep = Math.abs(stop - start) / Math.max(1, count)
  const power = Math.floor(Math.log10(rawStep || 1))
  const base = 10 ** power
  const error = rawStep / base
  const factor = error >= Math.sqrt(50) ? 10 : error >= Math.sqrt(10) ? 5 : error >= Math.sqrt(2) ? 2 : 1
  const step = factor * base
  const first = Math.ceil(start / step)
  const last = Math.floor(stop / step)
  return Array.from({ length: Math.max(0, last - first + 1) }, (_, index) => (first + index) * step)
}

function legendLayout(series: ChartSeries[], totalWidth: number) {
  const rows: { series: ChartSeries; width: number }[][] = []
  let row: { series: ChartSeries; width: number }[] = []
  let rowWidth = 0
  for (const item of series) {
    const width = Math.max(110, item.name.length * 8.6 + 46)
    if (row.length && rowWidth + 24 + width > totalWidth) {
      rows.push(row)
      row = []
      rowWidth = 0
    }
    row.push({ series: item, width })
    rowWidth += (row.length > 1 ? 24 : 0) + width
  }
  if (row.length) rows.push(row)
  return rows.flatMap((items, rowIndex) => {
    const width = items.reduce((sum, item, index) => sum + item.width + (index ? 24 : 0), 0)
    let cursor = (totalWidth - width) / 2
    return items.map((item, index) => {
      const entry = { series: item.series, x: cursor, y: rowIndex * 28 }
      cursor += item.width + (index < items.length - 1 ? 24 : 0)
      return entry
    })
  })
}

function formatFrequency(value: number) {
  if (Math.abs(value) >= 1000) return `${(value / 1000).toFixed(1)}k`
  if (Math.abs(value) >= 1) return value.toFixed(2).replace(/\.?0+$/, '')
  return value.toFixed(4).replace(/\.?0+$/, '')
}

function MetaPlotChart({ chart, svgRef }: { chart: ChartPayload; svgRef: React.RefObject<SVGSVGElement | null> }) {
  const [hover, setHover] = useState<ChartHover | null>(null)
  const plotRight = chartSize.width - chartMargin.right
  const plotBottom = chartSize.height - chartMargin.bottom
  const xDomain = chart.xDomain ?? [0, 1]
  const allValues = chart.series.flatMap((series) => series.values.flatMap((point) => [point.density, point.confidenceUp ?? point.density]))
  const yDomain = chart.yDomain ?? [-0.6, Math.max(1e-9, ...allValues.filter(Number.isFinite)) * 1.05]
  const xScale = (value: number) => chartMargin.left + (value - xDomain[0]) / (xDomain[1] - xDomain[0]) * (plotRight - chartMargin.left)
  const yScale = (value: number) => plotBottom - (value - yDomain[0]) / (yDomain[1] - yDomain[0]) * (plotBottom - chartMargin.top)
  const yTicks = tickValues(yDomain[0], yDomain[1])
  const legends = legendLayout(chart.series, plotRight - chartMargin.left)

  const showHover = (event: React.PointerEvent<SVGPathElement>, series: ChartSeries) => {
    const svg = event.currentTarget.ownerSVGElement
    const shell = svg?.closest('.rnaenrich-meta-chart-shell')
    if (!svg || !shell || !series.values.length) return
    const svgRect = svg.getBoundingClientRect()
    const shellRect = shell.getBoundingClientRect()
    const chartX = (event.clientX - svgRect.left) / svgRect.width * chartSize.width
    const dataX = xDomain[0] + (chartX - chartMargin.left) / (plotRight - chartMargin.left) * (xDomain[1] - xDomain[0])
    const point = series.values.reduce((nearest, candidate) => Math.abs(candidate.x - dataX) < Math.abs(nearest.x - dataX) ? candidate : nearest)
    const pointerX = event.clientX - shellRect.left
    const tooltipWidth = 210
    const left = pointerX + 24 + tooltipWidth > shellRect.width ? Math.max(12, pointerX - tooltipWidth - 22) : pointerX + 24
    const top = Math.max(12, Math.min(event.clientY - shellRect.top - 42, shellRect.height - 112))
    setHover({ series, point, left, top })
  }

  return <div className="rnaenrich-meta-chart-scroll">
    <div className="rnaenrich-meta-chart-shell">
      <svg ref={svgRef} className="rnaenrich-meta-chart" viewBox={`0 0 ${chartSize.width} ${chartSize.height}`} preserveAspectRatio="xMinYMin meet" role="img" aria-label={chart.title}>
        <rect width={chartSize.width} height={chartSize.height} fill="#fff" />
        <rect x={chartMargin.left} y={chartMargin.top} width={plotRight - chartMargin.left} height={plotBottom - chartMargin.top} fill="#ebebeb" />
        {[0, .25, .5, .75, 1].map((tick) => <line key={`x-${tick}`} x1={xScale(tick)} x2={xScale(tick)} y1={chartMargin.top} y2={plotBottom} stroke="rgba(255,255,255,.95)" strokeWidth="1.1" />)}
        {yTicks.map((tick) => <g key={`y-${tick}`}>
          <line x1={chartMargin.left} x2={plotRight} y1={yScale(tick)} y2={yScale(tick)} stroke="rgba(255,255,255,.95)" strokeWidth="1.1" />
          <line x1={chartMargin.left - 6} x2={chartMargin.left} y1={yScale(tick)} y2={yScale(tick)} stroke="#4c4c4c" strokeWidth="1.2" />
          <text x={chartMargin.left - 14} y={yScale(tick) + 5} textAnchor="end" fill="#4c4c4c" fontSize="14" fontFamily="sans-serif">{formatFrequency(tick)}</text>
        </g>)}
        <line x1={chartMargin.left} x2={chartMargin.left} y1={chartMargin.top} y2={plotBottom} stroke="#4c4c4c" strokeWidth="1.5" />
        <text x={chartMargin.left} y="50" fill="#111" fontSize="15" fontWeight="700" fontFamily="sans-serif">{chart.title}</text>
        <text transform="translate(30 250) rotate(-90)" fill="#111" fontSize="15" fontWeight="700" fontFamily="sans-serif">{chart.yLabel}</text>
        {chart.series.map((series) => <g key={series.name}>
          <path d={areaPath(series.values, xScale, yScale)} fill={series.color} opacity=".18" />
          {chart.showCI ? <>
            <path d={confidenceAreaPath(series.values, xScale, yScale)} fill={series.color} opacity=".12" />
            <path d={linePath(series.values, xScale, yScale, 'confidenceDown')} fill="none" stroke="#4f6db8" strokeWidth="1" opacity=".45" />
            <path d={linePath(series.values, xScale, yScale, 'confidenceUp')} fill="none" stroke="#111" strokeWidth="1" opacity=".45" />
          </> : null}
          <path d={linePath(series.values, xScale, yScale, 'density')} fill="none" stroke={series.color} strokeWidth="4" />
          <path data-export-ignore="true" d={linePath(series.values, xScale, yScale, 'density')} fill="none" stroke="transparent" strokeWidth="20" strokeLinecap="round" strokeLinejoin="round" className="rnaenrich-meta-hover-path" onPointerEnter={(event) => showHover(event, series)} onPointerMove={(event) => showHover(event, series)} onPointerLeave={() => setHover(null)} />
        </g>)}
        {chart.components.separators.map((separator) => <line key={separator.x} x1={xScale(separator.x)} x2={xScale(separator.x)} y1={chartMargin.top} y2={yScale(chart.components.trackBase ?? 0)} stroke="#111" strokeWidth="1.6" strokeDasharray="2 6" />)}
        {chart.components.segments.map((segment) => {
          const base = chart.components.trackBase ?? 0
          const minimum = ({ promoter: 2, tail: 2, utr5: 4, utr3: 4, cds: 6, ncrna: 6, rna: 6 } as Record<string, number>)[segment.component] ?? 3
          const height = Math.max(minimum, Math.abs(yScale(base - segment.height) - yScale(base + segment.height)))
          const center = yScale(base)
          const isFlank = segment.component === 'promoter' || segment.component === 'tail'
          return <g key={`${segment.component}-${segment.start}`}>
            <rect x={xScale(segment.start)} y={Math.round(center - height / 2) + .5} width={Math.max(1, xScale(segment.end) - xScale(segment.start))} height={Math.max(1, Math.round(height))} fill={isFlank ? '#111' : '#fff'} fillOpacity={isFlank ? 1 : segment.alpha} stroke="#111" strokeWidth="1.1" />
            {segment.label ? <text x={xScale(segment.mid)} y={yScale(chart.components.labelY ?? 0) + 12} textAnchor="middle" fill="#222" fontSize="14" fontWeight="700" fontFamily="sans-serif">{segment.label}</text> : null}
          </g>
        })}
        <g>
          {legends.map((entry) => <g transform={`translate(${chartMargin.left + entry.x} ${chartSize.height - 85 + entry.y})`} key={entry.series.name}>
            <rect x="0" y="-10" width="20" height="20" fill={entry.series.color} fillOpacity=".18" stroke={entry.series.color} strokeWidth="2" />
            <text x="34" y="0" dy=".35em" fill="#222" fontSize="15" fontWeight="700" fontFamily="sans-serif"><title>{entry.series.originalName}</title>{entry.series.name}</text>
          </g>)}
        </g>
      </svg>
      {hover ? <div className="rnaenrich-meta-tooltip" style={{ left: hover.left, top: hover.top }}>
        <strong>{hover.series.name}</strong>
        {hover.series.originalName !== hover.series.name ? <div><span>Source file</span><b>{hover.series.originalName}</b></div> : null}
        <div><span>Relative position</span><b>{(hover.point.x * 100).toFixed(1)}%</b></div>
        <div><span>{chart.yLabel || 'Frequency'}</span><b>{formatFrequency(hover.point.density)}</b></div>
      </div> : null}
    </div>
  </div>
}

function downloadBlob(filename: string, blob: Blob) {
  const url = URL.createObjectURL(blob)
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = filename
  document.body.appendChild(anchor)
  anchor.click()
  anchor.remove()
  window.setTimeout(() => URL.revokeObjectURL(url), 0)
}

function svgMarkup(svg: SVGSVGElement) {
  const clone = svg.cloneNode(true) as SVGSVGElement
  clone.querySelectorAll('[data-export-ignore="true"]').forEach((node) => node.remove())
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
  clone.setAttribute('width', String(chartSize.width))
  clone.setAttribute('height', String(chartSize.height))
  clone.setAttribute('font-family', 'sans-serif')
  clone.querySelectorAll('text, tspan').forEach((node) => node.setAttribute('font-family', 'sans-serif'))
  return clone
}

async function exportPng(svg: SVGSVGElement) {
  const source = new XMLSerializer().serializeToString(svgMarkup(svg))
  const url = URL.createObjectURL(new Blob([source], { type: 'image/svg+xml;charset=utf-8' }))
  const image = new Image()
  const blob = await new Promise<Blob>((resolve, reject) => {
    image.onerror = () => { URL.revokeObjectURL(url); reject(new Error('Unable to render the Meta Plot image.')) }
    image.onload = () => {
      const canvas = document.createElement('canvas')
      canvas.width = chartSize.width * 3
      canvas.height = chartSize.height * 3
      const context = canvas.getContext('2d')
      if (!context) { reject(new Error('Unable to create the export canvas.')); return }
      context.fillStyle = '#ffffff'
      context.fillRect(0, 0, canvas.width, canvas.height)
      context.drawImage(image, 0, 0, canvas.width, canvas.height)
      canvas.toBlob((value) => value ? resolve(value) : reject(new Error('Unable to encode PNG.')), 'image/png')
      URL.revokeObjectURL(url)
    }
    image.src = url
  })
  downloadBlob('meta_plot.png', blob)
}

async function exportPdf(svg: SVGSVGElement) {
  const pdf = new jsPDF({ orientation: 'landscape', unit: 'px', format: [chartSize.width, chartSize.height], compress: true })
  await svg2pdf(svgMarkup(svg), pdf, { x: 0, y: 0, width: chartSize.width, height: chartSize.height })
  downloadBlob('meta_plot.pdf', pdf.output('blob'))
}

function NumericField({ label, value, min, step = 1, onChange }: { label: string; value: number; min: number; step?: number; onChange: (value: number) => void }) {
  return <label className="rnaenrich-meta-field"><span>{label}</span><input type="number" min={min} step={step} value={value} onChange={(event) => onChange(Number(event.target.value))} /></label>
}

function StyledSelect({ label, value, options, disabled = false, onChange }: { label: string; value: string; options: SpeciesOption[]; disabled?: boolean; onChange: (value: string) => void }) {
  return <div className="rnaenrich-meta-field">
    <span>{label}</span>
    <Select.Root value={value} onValueChange={onChange} disabled={disabled}>
      <Select.Trigger className="rnaenrich-meta-select-trigger" aria-label={label}>
        <Select.Value />
        <Select.Icon className="rnaenrich-meta-select-icon"><ChevronDown size={15} /></Select.Icon>
      </Select.Trigger>
      <Select.Portal>
        <Select.Content className="rnaenrich-meta-select-content" position="popper" sideOffset={5} collisionPadding={12}>
          <Select.Viewport className="rnaenrich-meta-select-viewport">
            {options.map((option) => <Select.Item className="rnaenrich-meta-select-item" value={option.value} key={option.value}>
              <Select.ItemText>{option.label}</Select.ItemText>
              <Select.ItemIndicator className="rnaenrich-meta-select-indicator"><Check size={14} /></Select.ItemIndicator>
            </Select.Item>)}
          </Select.Viewport>
        </Select.Content>
      </Select.Portal>
    </Select.Root>
  </div>
}

export default function MetaPlot({ files }: { files: ParsedBed[] }) {
  const [speciesOptions, setSpeciesOptions] = useState<SpeciesOption[]>([])
  const [species, setSpecies] = useState('')
  const [controls, setControls] = useState<MetaPlotControls>(defaults)
  const [result, setResult] = useState<MetaPlotResult | null>(null)
  const [busy, setBusy] = useState(false)
  const [message, setMessage] = useState('Loading the annotation catalog…')
  const svgRef = useRef<SVGSVGElement | null>(null)

  useEffect(() => {
    let active = true
    fetch(`${apiBase}/species`).then(async (response) => {
      const body = await response.json() as { species?: SpeciesOption[]; error?: string }
      if (!response.ok) throw new Error(body.error || 'Unable to load species annotations.')
      if (!active) return
      const options = body.species ?? []
      setSpeciesOptions(options)
      setSpecies(options.find((option) => option.value === 'ara_TAIR10')?.value ?? options[0]?.value ?? '')
      setMessage(options.length ? '' : 'No species annotations are available.')
    }).catch((error: unknown) => { if (active) setMessage(error instanceof Error ? error.message : 'Unable to load the annotation catalog.') })
    return () => { active = false }
  }, [])

  useEffect(() => { setResult(null) }, [files])

  const speciesLabel = useMemo(() => speciesOptions.find((option) => option.value === species)?.label ?? species, [species, speciesOptions])
  const update = <K extends keyof MetaPlotControls>(key: K, value: MetaPlotControls[K]) => setControls((current) => ({ ...current, [key]: value }))

  const run = async () => {
    setMessage('')
    if (!files.length) { setMessage('Add at least one BED file before generating a profile.'); return }
    if (!species) { setMessage('Select a species annotation.'); return }
    setBusy(true)
    setResult(null)
    try {
      const response = await fetch(`${apiBase}/run`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          species,
          controls,
          files: files.map((file) => ({ name: file.name, content: file.records.map((record) => record.sourceLine).join('\n') })),
        }),
      })
      const body = await response.json() as MetaPlotResult & { error?: string }
      if (!response.ok) throw new Error(body.error || 'Meta Plot analysis failed.')
      setResult(body)
    } catch (error) {
      setMessage(error instanceof Error ? error.message : 'Meta Plot analysis failed.')
    } finally {
      setBusy(false)
    }
  }

  return <>
    <section className="rnaenrich-step rnaenrich-meta-settings">
      <div className="rnaenrich-step-heading"><span>02</span><div><h2>Configure Meta Plot</h2><p>Map genomic intervals to annotated transcript features and estimate their density across normalized transcript regions.</p></div></div>
      <div className="rnaenrich-meta-primary-grid">
        <StyledSelect label="Species annotation" value={species} options={speciesOptions} disabled={!speciesOptions.length} onChange={(value) => { setSpecies(value); setResult(null) }} />
        <StyledSelect label="Profile target" value={controls.pltTxType} options={[{ value: 'tx', label: 'Exon' }, { value: 'mrna', label: 'mRNA' }, { value: 'ncrna', label: 'lncRNA' }]} onChange={(value) => update('pltTxType', value as MetaPlotControls['pltTxType'])} />
        <StyledSelect label="Sampling scheme" value={controls.stSampleModle} options={[{ value: 'Equidistance', label: 'Equidistant' }, { value: 'random', label: 'Random' }]} onChange={(value) => update('stSampleModle', value as MetaPlotControls['stSampleModle'])} />
      </div>
      <div className="rnaenrich-meta-switches">
        {([
          ['headOrtail', 'Include promoter and tail'], ['txPrimaryOnly', 'Restrict to primary transcripts'],
          ['txlncrnaOverlapmrna', 'Retain overlapping lncRNA/mRNA'], ['mapFilterTranscript', 'Filter by source locus'],
          ['enableCI', 'Estimate confidence intervals'],
        ] as [keyof MetaPlotControls, string][]).map(([key, label]) => <label key={key}><input type="checkbox" checked={controls[key] === 'TRUE'} onChange={(event) => update(key, (event.target.checked ? 'TRUE' : 'FALSE') as never)} /><span>{label}</span></label>)}
      </div>
      <details className="rnaenrich-meta-advanced">
        <summary><SlidersHorizontal size={15} /> Advanced parameters</summary>
        <div className="rnaenrich-meta-numeric-grid">
          <NumericField label="Bases sampled per site" value={controls.stSampleNum} min={1} onChange={(value) => update('stSampleNum', value)} />
          <NumericField label="Maximum overlapping transcripts" value={controls.stAmblguity} min={0} onChange={(value) => update('stAmblguity', value)} />
          <NumericField label="Minimum 5′UTR length" value={controls.txfiveutrMinLength} min={0} onChange={(value) => update('txfiveutrMinLength', value)} />
          <NumericField label="Minimum CDS length" value={controls.txcdsMinLength} min={0} onChange={(value) => update('txcdsMinLength', value)} />
          <NumericField label="Minimum 3′UTR length" value={controls.txthreeutrMinLength} min={0} onChange={(value) => update('txthreeutrMinLength', value)} />
          <NumericField label="Minimum lncRNA length" value={controls.txlongNcrnaMinLength} min={0} onChange={(value) => update('txlongNcrnaMinLength', value)} />
          <NumericField label="Promoter length" value={controls.txpromoterLength} min={0} onChange={(value) => update('txpromoterLength', value)} />
          <NumericField label="Tail length" value={controls.txtailLength} min={0} onChange={(value) => update('txtailLength', value)} />
          <NumericField label="Density bandwidth adjustment" value={controls.adjust} min={0.01} step={0.1} onChange={(value) => update('adjust', value)} />
          <NumericField label="Transcript-overlap weight" value={controls.overlapIndex} min={0} step={0.1} onChange={(value) => update('overlapIndex', value)} />
          <NumericField label="Site-length weight" value={controls.siteLengthIndex} min={0} step={0.1} onChange={(value) => update('siteLengthIndex', value)} />
        </div>
      </details>
      <div className="rnaenrich-action-row">
        <div className="rnaenrich-action-summary"><strong>{files.length} BED file{files.length === 1 ? '' : 's'}</strong><span>{speciesLabel || 'No species selected'}</span></div>
        <button className="rnaenrich-button rnaenrich-button--primary" type="button" onClick={() => void run()} disabled={busy || !speciesOptions.length}>{busy ? 'Generating profile…' : 'Generate profile'}</button>
      </div>
      {message ? <div className="rnaenrich-message" role="status">{message}</div> : null}
    </section>

    {result ? <section className="rnaenrich-results" aria-label="Meta Plot results">
      <article className="rnaenrich-panel">
        <div className="rnaenrich-panel-heading rnaenrich-meta-result-heading"><div><LineChart size={18} /><h3>Transcript-scale distribution</h3></div><div className="rnaenrich-meta-exports">
          <button type="button" onClick={() => svgRef.current && void exportPng(svgRef.current)}><Download size={14} />PNG</button>
          <button type="button" onClick={() => svgRef.current && void exportPdf(svgRef.current)}><Download size={14} />PDF</button>
        </div></div>
        <div className="rnaenrich-meta-summary">
          <dl>
          <div><dt>Species</dt><dd>{speciesLabel}</dd></div>
          <div><dt>Input files</dt><dd>{result.summary.sample_count.toLocaleString()}</dd></div>
          <div><dt>Input intervals</dt><dd>{result.summary.interval_count.toLocaleString()}</dd></div>
          <div><dt>Normalized points</dt><dd>{result.summary.sampled_point_count.toLocaleString()}</dd></div>
          <div><dt>Mapped transcripts</dt><dd>{result.summary.transcript_hits.toLocaleString()}</dd></div>
          </dl>
        </div>
        <MetaPlotChart chart={result.chart} svgRef={svgRef} />
      </article>
    </section> : null}
  </>
}
