import { BarChart3, Download, FilePlus2, Grid3X3, Upload, X } from 'lucide-react'
import { jsPDF } from 'jspdf'
import { useEffect, useMemo, useRef, useState } from 'react'
import { svg2pdf } from 'svg2pdf.js'
import './RNAEnrich.css'
import MetaPlot from './MetaPlot'
import {
  compareBeds,
  exactMembership,
  exactSharedRecords,
  parseBed,
  recordKey,
  type BackgroundMode,
  type MatchMode,
  type PairResult,
  type ParsedBed,
} from './rnaenrich-core'

const appBase = import.meta.env.BASE_URL
const apiBase = `${appBase}RNAenrich/api`
const databaseCatalogUrl = `${appBase}RNAenrich/database-catalog.json`
const maxFiles = 6
const maxRecordsPerFile = 250_000

const demoBeds = Array.from({ length: maxFiles }, (_, fileIndex) => {
  const shared = Array.from({ length: 20 }, (_, index) => `chr1\t${100 + index * 100}\t${101 + index * 100}`)
  const partial = Array.from({ length: 8 }, (_, index) => (index + fileIndex) % 3 === 0 ? `chr2\t${500 + index * 120}\t${501 + index * 120}` : null).filter((line): line is string => line !== null)
  const unique = Array.from({ length: 12 }, (_, index) => `chr${fileIndex + 3}\t${1000 + index * 80}\t${1001 + index * 80}`)
  return {
    name: `rnaenrich_demo_${fileIndex + 1}.bed`,
    text: [...shared, ...partial, ...unique].join('\n'),
  }
})

type Metric = 'jaccard' | 'overlap' | 'fdr'
type AnalysisTool = 'overlap' | 'meta'
type DatabaseFileOption = { key: string; label: string; accession?: string; sample?: string; treatment?: string; caller?: string; site_count?: number; file_url?: string }
type DatabaseSourceOption = { id: string; label: string; species: string; species_label: string; file_count: number; files: DatabaseFileOption[] }
type DatabaseSelectOption = { value: string; label: string }

function DatabaseSelect({ value, options, placeholder, disabled, onChange }: { value: string; options: DatabaseSelectOption[]; placeholder: string; disabled?: boolean; onChange: (value: string) => void }) {
  const [open, setOpen] = useState(false)
  const rootRef = useRef<HTMLDivElement>(null)
  const selected = options.find((option) => option.value === value)
  useEffect(() => {
    if (!open) return
    const closeOnOutsideClick = (event: MouseEvent) => {
      if (!rootRef.current?.contains(event.target as Node)) setOpen(false)
    }
    const closeOnEscape = (event: KeyboardEvent) => { if (event.key === 'Escape') setOpen(false) }
    document.addEventListener('mousedown', closeOnOutsideClick)
    document.addEventListener('keydown', closeOnEscape)
    return () => { document.removeEventListener('mousedown', closeOnOutsideClick); document.removeEventListener('keydown', closeOnEscape) }
  }, [open])
  return <div className="rnaenrich-custom-select" ref={rootRef}>
    <button className="rnaenrich-custom-select-trigger" type="button" disabled={disabled} aria-haspopup="listbox" aria-expanded={open} onClick={() => setOpen((current) => !current)}>
      <span>{selected?.label ?? placeholder}</span><i aria-hidden="true" />
    </button>
    {open ? <div className="rnaenrich-custom-select-menu" role="listbox">
      {options.length ? options.map((option) => <button className={option.value === value ? 'is-selected' : ''} type="button" role="option" aria-selected={option.value === value} key={option.value} onClick={() => { onChange(option.value); setOpen(false) }}>{option.label}</button>) : <span className="rnaenrich-custom-select-empty">No options available</span>}
    </div> : null}
  </div>
}

const formatPercent = (value: number) => `${(value * 100).toFixed(1)}%`
const formatProbability = (value: number | null) => {
  if (value === null) return '—'
  if (value < 0.001) return value.toExponential(2)
  return value.toFixed(3)
}

const safeFilename = (value: string) => value.replace(/\.[^.]+$/, '').replace(/[^a-zA-Z0-9._-]+/g, '_')
const compactFilename = (value: string, index: number) => {
  const compact = value.length > 24 ? `${value.slice(0, 11)}…${value.slice(-10)}` : value
  return `${Math.max(0, index) + 1} · ${compact}`
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

function downloadText(filename: string, content: string, type = 'text/plain;charset=utf-8') {
  downloadBlob(filename, new Blob([content], { type }))
}

const delimitedCell = (value: string | number | null, delimiter: ',' | '\t') => {
  const text = value === null ? '' : String(value)
  return text.includes('"') || text.includes('\n') || text.includes('\r') || text.includes(delimiter) ? `"${text.replace(/"/g, '""')}"` : text
}

function downloadDelimited(filename: string, rows: (string | number | null)[][], delimiter: ',' | '\t') {
  downloadText(filename, `${rows.map((row) => row.map((value) => delimitedCell(value, delimiter)).join(delimiter)).join('\n')}\n`, delimiter === ',' ? 'text/csv;charset=utf-8' : 'text/plain;charset=utf-8')
}

type SvgExport = { markup: string; width: number; height: number }

const escapeXml = (value: string | number) => String(value).replace(/[&<>"']/g, (character) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&apos;' })[character] ?? character)
const exportLabel = (value: string, maximum = 30) => value.length > maximum ? `${value.slice(0, maximum - 1)}…` : value

function mixHex(background: string, foreground: string, amount: number) {
  const channel = (color: string, offset: number) => Number.parseInt(color.slice(offset, offset + 2), 16)
  const mixed = [1, 3, 5].map((offset) => Math.round(channel(background, offset) * (1 - amount) + channel(foreground, offset) * amount).toString(16).padStart(2, '0'))
  return `#${mixed.join('')}`
}

function matrixValue(results: PairResult[], metric: Metric, row: ParsedBed, column: ParsedBed) {
  if (row.name === column.name) return metric === 'jaccard' ? 1 : metric === 'overlap' ? row.records.length : null
  const pair = pairFor(results, row.name, column.name)
  if (!pair) return null
  if (metric === 'jaccard') return pair.jaccard
  if (metric === 'overlap') return pair.overlapCount
  return pair.fdr === null ? null : -Math.log10(Math.max(pair.fdr, 1e-300))
}

function buildMatrixSvg(files: ParsedBed[], results: PairResult[], metric: Metric): SvgExport {
  const width = Math.max(900, 260 + files.length * 160)
  const left = 250
  const top = 100
  const cellWidth = (width - left - 30) / files.length
  const cellHeight = 60
  const height = top + files.length * cellHeight + 35
  const label = metric === 'jaccard' ? 'Jaccard' : metric === 'overlap' ? 'Shared sites' : '−log10 FDR'
  const values = files.flatMap((row) => files.map((column) => matrixValue(results, metric, row, column))).filter((value): value is number => value !== null)
  const maximum = metric === 'jaccard' ? 1 : Math.max(1, ...values)
  const columns = files.map((file, index) => `<text x="${left + (index + .5) * cellWidth}" y="78" text-anchor="middle" fill="#607778" font-size="12" font-weight="700">${escapeXml(exportLabel(file.name))}</text>`).join('')
  const rows = files.map((row, rowIndex) => {
    const rowLabel = `<text x="20" y="${top + (rowIndex + .5) * cellHeight + 4}" fill="#607778" font-size="12" font-weight="700">${escapeXml(exportLabel(row.name))}</text>`
    const cells = files.map((column, columnIndex) => {
      const value = matrixValue(results, metric, row, column)
      const intensity = value === null ? 0 : Math.max(.08, value / maximum)
      const fill = value === null ? '#f0f3f1' : mixHex('#edf4f1', '#397b73', intensity * .88)
      const text = value === null ? '—' : metric === 'jaccard' ? value.toFixed(3) : metric === 'fdr' ? value.toFixed(2) : value.toLocaleString()
      return `<rect x="${left + columnIndex * cellWidth}" y="${top + rowIndex * cellHeight}" width="${cellWidth}" height="${cellHeight}" fill="${fill}" stroke="#ffffff" stroke-width="2"/><text x="${left + (columnIndex + .5) * cellWidth}" y="${top + (rowIndex + .5) * cellHeight + 5}" text-anchor="middle" fill="#123233" font-size="13" font-weight="800">${escapeXml(text)}</text>`
    }).join('')
    return rowLabel + cells
  }).join('')
  return { width, height, markup: `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="sans-serif"><rect width="100%" height="100%" fill="#ffffff"/><text x="20" y="32" fill="#143536" font-size="18" font-weight="800">Pairwise similarity matrix · ${escapeXml(label)}</text>${columns}${rows}</svg>` }
}

function buildUpSetSvg(files: ParsedBed[], combinations: ReturnType<typeof exactMembership>): SvgExport {
  const width = Math.max(900, 230 + combinations.length * 90)
  const left = 220
  const barBottom = 210
  const barMaximumHeight = 120
  const rowTop = 245
  const rowHeight = 34
  const height = rowTop + files.length * rowHeight + 30
  const plotWidth = width - left - 30
  const columnWidth = plotWidth / Math.max(1, combinations.length)
  const maximum = Math.max(1, ...combinations.map((item) => item.count))
  const fileRows = files.map((file, index) => `<line x1="20" x2="${width - 30}" y1="${rowTop + (index + 1) * rowHeight}" y2="${rowTop + (index + 1) * rowHeight}" stroke="#edf1ef"/><text x="20" y="${rowTop + index * rowHeight + 22}" fill="#607778" font-size="12" font-weight="700">${escapeXml(exportLabel(file.name))}</text>`).join('')
  const columns = combinations.map((combination, index) => {
    const x = left + (index + .5) * columnWidth
    const barHeight = Math.max(3, combination.count / maximum * barMaximumHeight)
    const activeIndexes = files.flatMap((file, fileIndex) => combination.members.includes(file.name) ? [fileIndex] : [])
    const first = Math.min(...activeIndexes)
    const last = Math.max(...activeIndexes)
    const line = last > first ? `<line x1="${x}" x2="${x}" y1="${rowTop + first * rowHeight + rowHeight / 2}" y2="${rowTop + last * rowHeight + rowHeight / 2}" stroke="#397b73" stroke-width="2"/>` : ''
    const dots = files.map((file, fileIndex) => `<circle cx="${x}" cy="${rowTop + fileIndex * rowHeight + rowHeight / 2}" r="${combination.members.includes(file.name) ? 6 : 4}" fill="${combination.members.includes(file.name) ? '#397b73' : '#d8e1de'}"/>`).join('')
    return `<text x="${x}" y="${barBottom - barHeight - 10}" text-anchor="middle" fill="#143536" font-size="12" font-weight="800">${combination.count.toLocaleString()}</text><rect x="${x - 12}" y="${barBottom - barHeight}" width="24" height="${barHeight}" rx="3" fill="#c66f62"/>${line}${dots}`
  }).join('')
  return { width, height, markup: `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="sans-serif"><rect width="100%" height="100%" fill="#ffffff"/><text x="20" y="32" fill="#143536" font-size="18" font-weight="800">Set intersections</text><text x="20" y="54" fill="#607778" font-size="12" font-weight="700">Count</text><line x1="20" x2="${width - 30}" y1="${barBottom}" y2="${barBottom}" stroke="#d7e2dd"/>${fileRows}${columns}</svg>` }
}

async function svgToPng(plot: SvgExport, filename: string) {
  const url = URL.createObjectURL(new Blob([plot.markup], { type: 'image/svg+xml;charset=utf-8' }))
  const image = new Image()
  const blob = await new Promise<Blob>((resolve, reject) => {
    image.onerror = () => { URL.revokeObjectURL(url); reject(new Error('Unable to render the exported image.')) }
    image.onload = () => {
      const scale = 2
      const canvas = document.createElement('canvas')
      canvas.width = plot.width * scale
      canvas.height = plot.height * scale
      const context = canvas.getContext('2d')
      if (!context) { URL.revokeObjectURL(url); reject(new Error('Unable to create the export canvas.')); return }
      context.scale(scale, scale)
      context.fillStyle = '#ffffff'
      context.fillRect(0, 0, plot.width, plot.height)
      context.drawImage(image, 0, 0, plot.width, plot.height)
      canvas.toBlob((value) => value ? resolve(value) : reject(new Error('Unable to encode the exported image.')), 'image/png')
      URL.revokeObjectURL(url)
    }
    image.src = url
  })
  downloadBlob(filename, blob)
}

async function svgToPdf(plot: SvgExport, filename: string) {
  const documentNode = new DOMParser().parseFromString(plot.markup, 'image/svg+xml')
  if (documentNode.querySelector('parsererror')) throw new Error('Unable to parse the exported figure.')
  const svg = documentNode.documentElement as unknown as SVGSVGElement
  const pdf = new jsPDF({ orientation: plot.width >= plot.height ? 'landscape' : 'portrait', unit: 'px', format: [plot.width, plot.height], compress: true })
  await svg2pdf(svg, pdf, { x: 0, y: 0, width: plot.width, height: plot.height })
  downloadBlob(filename, pdf.output('blob'))
}

function pairFor(results: PairResult[], left: string, right: string) {
  return results.find((pair) => (pair.left === left && pair.right === right) || (pair.left === right && pair.right === left))
}

function Heatmap({ files, results, metric }: { files: ParsedBed[]; results: PairResult[]; metric: Metric }) {
  const cellValue = (row: ParsedBed, column: ParsedBed) => {
    if (row.name === column.name) return metric === 'jaccard' ? 1 : metric === 'overlap' ? row.records.length : null
    const pair = pairFor(results, row.name, column.name)
    if (!pair) return null
    if (metric === 'jaccard') return pair.jaccard
    if (metric === 'overlap') return pair.overlapCount
    return pair.fdr === null ? null : -Math.log10(Math.max(pair.fdr, 1e-300))
  }
  const values = files.flatMap((row) => files.map((column) => cellValue(row, column))).filter((value): value is number => value !== null)
  const maxValue = metric === 'jaccard' ? 1 : Math.max(1, ...values)

  return <div className="rnaenrich-heatmap-scroll">
    <div className="rnaenrich-heatmap" style={{ gridTemplateColumns: `minmax(7rem, 1.25fr) repeat(${files.length}, minmax(4.25rem, 1fr))` }}>
      <div className="rnaenrich-heatmap-corner">File</div>
      {files.map((file, index) => <div className="rnaenrich-heatmap-axis rnaenrich-heatmap-axis--top" title={file.name} key={file.name}>{compactFilename(file.name, index)}</div>)}
      {files.map((row, rowIndex) => <div className="rnaenrich-heatmap-row" key={row.name}>
        <div className="rnaenrich-heatmap-axis" title={row.name}>{compactFilename(row.name, rowIndex)}</div>
        {files.map((column) => {
          const value = cellValue(row, column)
          const intensity = value === null ? 0 : Math.max(0.08, value / maxValue)
          const label = value === null ? '—' : metric === 'jaccard' ? value.toFixed(3) : metric === 'fdr' ? value.toFixed(2) : value.toLocaleString()
          return <div
            className={`rnaenrich-heatmap-cell${value === null ? ' is-empty' : ''}`}
            style={value === null ? undefined : { backgroundColor: `color-mix(in srgb, #397b73 ${Math.round(intensity * 88)}%, #edf4f1)` }}
            title={`${row.name} × ${column.name}: ${label}`}
            key={column.name}
          >{label}</div>
        })}
      </div>)}
    </div>
  </div>
}

function UpSetPlot({ files, combinations }: { files: ParsedBed[]; combinations: ReturnType<typeof exactMembership> }) {
  const maximum = Math.max(1, ...combinations.map((item) => item.count))
  return <div className="rnaenrich-upset-scroll">
    <div className="rnaenrich-upset-plot" style={{ gridTemplateColumns: `minmax(10rem, 13rem) repeat(${combinations.length}, minmax(3.8rem, 1fr))` }}>
      <div className="rnaenrich-upset-labels" style={{ gridTemplateRows: `1.5rem 7rem repeat(${files.length}, 1.7rem)` }}>
        <strong>Count</strong>
        <div aria-hidden="true" />
        {files.map((file, index) => <span title={file.name} key={file.name}>{compactFilename(file.name, index)}</span>)}
      </div>
      {combinations.map((combination) => {
        const activeIndexes = files.flatMap((file, index) => combination.members.includes(file.name) ? [index] : [])
        const first = Math.min(...activeIndexes)
        const last = Math.max(...activeIndexes)
        return <div className="rnaenrich-upset-column" style={{ gridTemplateRows: `1.5rem 7rem repeat(${files.length}, 1.7rem)` }} title={`${combination.label}: ${combination.count.toLocaleString()}`} key={combination.label}>
          <strong>{combination.count.toLocaleString()}</strong>
          <div className="rnaenrich-upset-bar"><i style={{ height: `${Math.max(3, combination.count / maximum * 100)}%` }} /></div>
          <div className="rnaenrich-upset-dots" style={{ gridTemplateRows: `repeat(${files.length}, 1.7rem)` }}>
            {last > first ? <i className="rnaenrich-upset-line" style={{ top: `${first * 1.7 + .75}rem`, height: `${(last - first) * 1.7}rem` }} /> : null}
            {files.map((file) => <span className={combination.members.includes(file.name) ? 'is-active' : ''} key={file.name} />)}
          </div>
        </div>
      })}
    </div>
  </div>
}

export default function RNAEnrich() {
  const [analysisTool, setAnalysisTool] = useState<AnalysisTool>('overlap')
  const [files, setFiles] = useState<ParsedBed[]>([])
  const [mode, setMode] = useState<MatchMode>('exact')
  const [backgroundMode, setBackgroundMode] = useState<BackgroundMode>('none')
  const [background, setBackground] = useState<ParsedBed | null>(null)
  const [results, setResults] = useState<PairResult[] | null>(null)
  const [metric, setMetric] = useState<Metric>('jaccard')
  const [message, setMessage] = useState('')
  const [databaseCatalog, setDatabaseCatalog] = useState<DatabaseSourceOption[]>([])
  const [databaseSpecies, setDatabaseSpecies] = useState('')
  const [databaseSource, setDatabaseSource] = useState('')
  const [databaseFiles, setDatabaseFiles] = useState<string[]>([])
  const [databaseEnabled, setDatabaseEnabled] = useState(false)
  const [databaseLoading, setDatabaseLoading] = useState(false)
  const [databaseAdding, setDatabaseAdding] = useState(false)
  const [uploadedFileNames, setUploadedFileNames] = useState<string[]>([])

  useEffect(() => {
    if (!databaseEnabled || databaseCatalog.length) return
    let active = true
    setDatabaseLoading(true)
    fetch(databaseCatalogUrl).then(async (response) => {
      const body = await response.json() as { databases?: DatabaseSourceOption[]; error?: string }
      if (!response.ok) throw new Error(body.error || 'Unable to load the database catalog.')
      if (!active) return
      const databases = body.databases ?? []
      setDatabaseCatalog(databases)
      setDatabaseSpecies('')
      setDatabaseSource('')
    }).catch((error: unknown) => {
      if (active) setMessage(error instanceof Error ? error.message : 'Unable to load the database catalog.')
    }).finally(() => { if (active) setDatabaseLoading(false) })
    return () => { active = false }
  }, [databaseCatalog.length, databaseEnabled])

  const allCombinations = useMemo(() => mode === 'exact' && results ? exactMembership(files) : [], [files, mode, results])
  const combinations = allCombinations.slice(0, 20)
  const unionSize = useMemo(() => new Set(files.flatMap((file) => file.records.map(recordKey))).size, [files])
  const databaseSpeciesOptions = useMemo(() => [...new Map(databaseCatalog.map((entry) => [entry.species, entry.species_label])).entries()], [databaseCatalog])
  const databaseSourceOptions = useMemo(() => databaseCatalog.filter((entry) => entry.species === databaseSpecies), [databaseCatalog, databaseSpecies])
  const selectedDatabaseSource = databaseSourceOptions.find((entry) => entry.id === databaseSource) ?? null

  const closeDatabaseModal = () => {
    setDatabaseEnabled(false)
    setDatabaseSpecies('')
    setDatabaseSource('')
    setDatabaseFiles([])
  }

  const selectDatabaseSpecies = (value: string) => {
    setDatabaseSpecies(value)
    setDatabaseSource('')
    setDatabaseFiles([])
  }

  const selectDatabaseSource = (value: string) => {
    setDatabaseSource(value)
    setDatabaseFiles([])
  }

  const addDatabaseFiles = async () => {
    if (databaseAdding) return
    if (!selectedDatabaseSource || !databaseFiles.length) {
      setMessage('Select a database and at least one BED file.')
      return
    }
    const available = maxFiles - files.length
    if (available <= 0) {
      setMessage(`RNAenrich accepts up to ${maxFiles} files per run.`)
      return
    }
    setDatabaseAdding(true)
    setMessage('Adding selected database BED files…')
    const keys = databaseFiles.slice(0, available)
    const existingNames = new Set(files.map((file) => file.name))
    try {
      const parsed = await Promise.all(keys.map(async (key) => {
        const fileOption = selectedDatabaseSource.files.find((file) => file.key === key)
        const response = fileOption?.file_url
          ? await fetch(`${appBase}${fileOption.file_url}`)
          : await fetch(`${apiBase}/database/file`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ species: selectedDatabaseSource.species, source_id: selectedDatabaseSource.id, file_key: key }),
          })
        if (!response.ok) throw new Error(`Unable to load ${key}.`)
        const originalName = fileOption?.label || key
        let name = originalName
        let suffix = 2
        while (existingNames.has(name)) {
          const stem = originalName.replace(/(\.[^.]+)?$/, '')
          const extension = originalName.match(/\.[^.]+$/)?.[0] ?? ''
          name = `${stem}-${suffix}${extension}`
          suffix += 1
        }
        existingNames.add(name)
        return parseBed(name, await response.text())
      }))
      setFiles((current) => [...current, ...parsed])
      setUploadedFileNames((current) => [...current, ...parsed.map((file) => file.name)])
      setDatabaseFiles([])
      setResults(null)
      closeDatabaseModal()
      setMessage('')
    } catch (error) {
      setMessage(error instanceof Error ? error.message : 'Unable to load database BED files.')
    } finally {
      setDatabaseAdding(false)
    }
  }

  const matrixRows = () => {
    const label = metric === 'jaccard' ? 'Jaccard' : metric === 'overlap' ? 'Shared sites' : '-log10 FDR'
    const value = (row: ParsedBed, column: ParsedBed) => {
      if (row.name === column.name) return metric === 'jaccard' ? 1 : metric === 'overlap' ? row.records.length : null
      const pair = pairFor(results ?? [], row.name, column.name)
      if (!pair) return null
      if (metric === 'jaccard') return pair.jaccard
      if (metric === 'overlap') return pair.overlapCount
      return pair.fdr === null ? null : -Math.log10(Math.max(pair.fdr, 1e-300))
    }
    return [[label, ...files.map((file) => file.name)], ...files.map((row) => [row.name, ...files.map((column) => value(row, column))])]
  }

  const intersectionRows = () => [
    ['Combination', 'Member files', 'Count'],
    ...allCombinations.map((combination) => [combination.label, combination.members.join('; '), combination.count]),
  ]

  const pairRows = () => [
    ['First file', 'Second file', 'First sites', 'Second sites', 'Shared sites', 'First matched', 'Second matched', 'Jaccard', 'P-value', 'FDR'],
    ...(results ?? []).map((pair) => [pair.left, pair.right, pair.leftSize, pair.rightSize, pair.overlapCount, pair.leftRate, pair.rightRate, pair.jaccard, pair.pValue, pair.fdr]),
  ]

  const exportMatrixCsv = () => downloadDelimited(`rnaenrich_similarity_${metric}.csv`, matrixRows(), ',')
  const exportMatrixTxt = () => downloadDelimited(`rnaenrich_similarity_${metric}.txt`, matrixRows(), '\t')
  const exportIntersectionsCsv = () => downloadDelimited('rnaenrich_set_intersections.csv', intersectionRows(), ',')
  const exportIntersectionsTxt = () => downloadDelimited('rnaenrich_set_intersections.txt', intersectionRows(), '\t')
  const exportPairsCsv = () => downloadDelimited('rnaenrich_pairwise_metrics.csv', pairRows(), ',')
  const exportPairsTxt = () => downloadDelimited('rnaenrich_pairwise_metrics.txt', pairRows(), '\t')

  const loadDemo = () => {
    setFiles(demoBeds.map(({ name, text }) => parseBed(name, text)))
    setUploadedFileNames([])
    setBackground(null)
    setBackgroundMode('none')
    setMode('exact')
    setResults(null)
    setMetric('jaccard')
    setMessage('')
  }

  const exportFigure = (plot: SvgExport, filename: string, format: 'png' | 'pdf') => {
    setMessage('')
    const task = format === 'png' ? svgToPng(plot, `${filename}.png`) : svgToPdf(plot, `${filename}.pdf`)
    void task.catch((error: unknown) => setMessage(error instanceof Error ? error.message : `Unable to export ${format.toUpperCase()}.`))
  }

  const addFiles = async (selected: FileList | null) => {
    if (!selected?.length) return
    setMessage('')
    const available = maxFiles - files.length
    const incoming = [...selected].slice(0, available)
    const existingNames = new Set(files.map((file) => file.name))
    const parsed = await Promise.all(incoming.map(async (file) => {
      let name = file.name
      let suffix = 2
      while (existingNames.has(name)) {
        const stem = file.name.replace(/(\.[^.]+)?$/, '')
        const extension = file.name.match(/\.[^.]+$/)?.[0] ?? ''
        name = `${stem}-${suffix}${extension}`
        suffix += 1
      }
      existingNames.add(name)
      return parseBed(name, await file.text())
    }))
    const oversized = parsed.find((file) => file.records.length > maxRecordsPerFile)
    if (oversized) {
      setMessage(`${oversized.name} contains more than ${maxRecordsPerFile.toLocaleString()} unique intervals. Split or filter it before analysis.`)
      return
    }
    setFiles((current) => [...current, ...parsed])
    setUploadedFileNames((current) => [...current, ...parsed.map((file) => file.name)])
    setResults(null)
    if (selected.length > available) setMessage(`Only the first ${available} file(s) were added; RNAenrich accepts up to ${maxFiles} files per run.`)
  }

  const uploadBackground = async (selected: File | undefined) => {
    if (!selected) return
    const parsed = parseBed(selected.name, await selected.text())
    if (parsed.records.length > maxRecordsPerFile) {
      setMessage(`Background contains more than ${maxRecordsPerFile.toLocaleString()} unique intervals.`)
      return
    }
    setBackground(parsed)
    setResults(null)
    setMessage('')
  }

  const runAnalysis = () => {
    setMessage('')
    if (files.length < 2) {
      setMessage('Add at least two BED files before running overlap analysis.')
      return
    }
    const empty = files.find((file) => file.records.length === 0)
    if (empty) {
      setMessage(`${empty.name} has no valid BED records.`)
      return
    }
    let backgroundSize: number | null = null
    if (mode === 'exact' && backgroundMode === 'union') backgroundSize = unionSize
    if (mode === 'exact' && backgroundMode === 'upload') {
      if (!background?.records.length) {
        setMessage('Upload a valid background BED or choose another background option.')
        return
      }
      const backgroundKeys = new Set(background.records.map(recordKey))
      const outside = files.flatMap((file) => file.records.filter((record) => !backgroundKeys.has(recordKey(record))).slice(0, 1).map((record) => `${file.name}: ${recordKey(record)}`))
      if (outside.length) {
        setMessage(`The background does not contain every input site. First unmatched site: ${outside[0]}`)
        return
      }
      backgroundSize = background.records.length
    }
    setResults(compareBeds(files, mode, backgroundSize))
    setMetric('jaccard')
  }

  return <div className="rnaenrich-page">
    <header className="rnaenrich-header">
      <a className="rnaenrich-brand" href={appBase}>RNA<span>meta</span></a>
      <nav aria-label="RNAenrich navigation">
        <a href={appBase}>Home</a>
        <a href={`${appBase}#news`}>News</a>
        <a href="https://github.com/Rwtkc/RNAmeta_R" target="_blank" rel="noreferrer">GitHub</a>
      </nav>
    </header>

    <main className="rnaenrich-main">
      <section className="rnaenrich-intro" aria-label="About RNAenrich">
        <div className="rnaenrich-intro-copy">
          <strong>RNAenrich evaluates concordance among genomic interval sets and characterizes their distributions across normalized transcript regions.</strong>
          <p>Use <strong>Overlap analysis</strong> to quantify exact-coordinate matches or interval intersections across multiple BED files. Use <strong>Meta Plot</strong> to map genomic intervals to annotated exon, mRNA, or lncRNA structures and compare their transcript-scale density profiles.</p>
        </div>
      </section>

      <div className="rnaenrich-tool-tabs" role="tablist" aria-label="RNAenrich analysis">
        <button type="button" role="tab" aria-selected={analysisTool === 'overlap'} className={analysisTool === 'overlap' ? 'is-active' : ''} onClick={() => setAnalysisTool('overlap')}>Overlap analysis</button>
        <button type="button" role="tab" aria-selected={analysisTool === 'meta'} className={analysisTool === 'meta' ? 'is-active' : ''} onClick={() => setAnalysisTool('meta')}>Meta Plot</button>
      </div>

      <section className="rnaenrich-workbench" aria-label="RNAenrich analysis settings">
        <div className="rnaenrich-step">
          <div className="rnaenrich-step-heading"><span>01</span><div><h2>Add BED files</h2><p>Upload {analysisTool === 'overlap' ? '2' : '1'}–{maxFiles} BED files. Genomic coordinates are read from the first three columns; invalid records are excluded and duplicate coordinates are collapsed.</p></div></div>
          <label className="rnaenrich-empty-upload">
            <Upload size={20} aria-hidden="true" />
            <strong>{files.length ? 'Add BED files' : 'Choose BED files'}</strong>
            <input type="file" multiple accept=".bed,.txt,text/plain" onChange={(event) => { void addFiles(event.target.files); event.target.value = '' }} />
          </label>
          <div className="rnaenrich-database-picker">
            <div className="rnaenrich-database-heading">
              <div><strong>Select from RNAmeta databases</strong><span>Add curated BED sites after uploading your own BED file.</span></div>
              <button className="rnaenrich-database-toggle" type="button" aria-expanded={databaseEnabled} disabled={!uploadedFileNames.length} onClick={() => setDatabaseEnabled(true)}>Use database BED files</button>
            </div>
          </div>
          {databaseEnabled ? <div className="rnaenrich-modal-backdrop" role="presentation" onMouseDown={(event) => { if (!databaseAdding && event.target === event.currentTarget) closeDatabaseModal() }}>
            <section className="rnaenrich-database-modal" role="dialog" aria-modal="true" aria-labelledby="rnaenrich-database-title">
              <div className="rnaenrich-database-modal-heading"><div><strong id="rnaenrich-database-title">Select RNAmeta database sites</strong><span>Choose a species, database, and BED files to add.</span></div><button className="rnaenrich-modal-close" type="button" aria-label="Close database selection" onClick={closeDatabaseModal} disabled={databaseAdding}><X size={18} /></button></div>
              <div className="rnaenrich-database-fields">
                <label><span>Species</span><DatabaseSelect value={databaseSpecies} disabled={databaseLoading || !databaseSpeciesOptions.length} placeholder="Select species" options={databaseSpeciesOptions.map(([value, label]) => ({ value, label }))} onChange={selectDatabaseSpecies} /></label>
                <label><span>Database</span><DatabaseSelect value={databaseSource} disabled={databaseLoading || !databaseSourceOptions.length} placeholder="Select database" options={databaseSourceOptions.map((entry) => ({ value: entry.id, label: `${entry.label} · ${entry.file_count} files` }))} onChange={selectDatabaseSource} /></label>
              </div>
              <div className="rnaenrich-database-files"><span>Database BED files</span><div className="rnaenrich-database-table-wrap"><table className="rnaenrich-database-table"><thead><tr><th aria-label="Select" /><th>BED file</th><th>Sites</th><th>Experiment</th></tr></thead><tbody>{selectedDatabaseSource?.files.map((file) => { const checked = databaseFiles.includes(file.key); const experiment = [file.sample, file.treatment, file.caller].filter(Boolean).join(' · ') || '—'; return <tr className={checked ? 'is-selected' : ''} key={file.key} onClick={() => setDatabaseFiles((current) => current.includes(file.key) ? current.filter((key) => key !== file.key) : [...current, file.key])}><td><input type="checkbox" checked={checked} aria-label={`Select ${file.label}`} onChange={() => undefined} /></td><td title={file.label}>{file.label}</td><td>{typeof file.site_count === 'number' ? file.site_count.toLocaleString() : '—'}</td><td title={experiment}>{experiment}</td></tr> })}</tbody></table></div></div>
              <div className="rnaenrich-database-footer"><small>{databaseAdding ? 'Reading and parsing selected files…' : databaseLoading ? 'Loading database catalog…' : selectedDatabaseSource ? `${selectedDatabaseSource.files.length} BED files available.` : 'Select a database.'}</small><div><button className="rnaenrich-button rnaenrich-button--quiet" type="button" onClick={closeDatabaseModal} disabled={databaseAdding}>Cancel</button><button className="rnaenrich-button rnaenrich-button--primary" type="button" onClick={() => void addDatabaseFiles()} disabled={databaseAdding || !selectedDatabaseSource || !databaseFiles.length}>{databaseAdding ? 'Adding…' : 'Add selected BED files'}</button></div></div>
            </section>
          </div> : null}
          {files.length ? <div className="rnaenrich-file-list">
            {files.map((file) => <article key={file.name}>
              <FilePlus2 size={17} aria-hidden="true" />
              <div><strong title={file.name}>{file.name}</strong><small>{file.records.length.toLocaleString()} unique · {file.duplicateCount} duplicate · {file.invalidLines.length} invalid</small></div>
              <button type="button" aria-label={`Remove ${file.name}`} onClick={() => { setFiles((current) => current.filter((item) => item.name !== file.name)); setUploadedFileNames((current) => current.filter((name) => name !== file.name)); setResults(null) }}><X size={16} /></button>
            </article>)}
          </div> : null}
        </div>

        {analysisTool === 'overlap' ? <div className="rnaenrich-step">
          <div className="rnaenrich-step-heading"><span>02</span><div><h2>Define the comparison</h2><p>Quantify concordance using identical genomic coordinates or intervals that intersect by at least one base pair.</p></div></div>
          <fieldset className="rnaenrich-choice-grid">
            <legend>Match mode</legend>
            <label className={mode === 'exact' ? 'is-selected' : ''}><input type="radio" name="mode" checked={mode === 'exact'} onChange={() => { setMode('exact'); setResults(null) }} /><strong>Exact match</strong><small>Chromosome, start, and end must be identical.</small></label>
            <label className={mode === 'interval' ? 'is-selected' : ''}><input type="radio" name="mode" checked={mode === 'interval'} onChange={() => { setMode('interval'); setResults(null) }} /><strong>Interval overlap</strong><small>Intervals on the same chromosome must share at least 1 bp.</small></label>
          </fieldset>
          <fieldset className="rnaenrich-choice-grid rnaenrich-choice-grid--three" disabled={mode === 'interval'}>
            <legend>Statistical background {mode === 'interval' ? '— not used for interval mode' : ''}</legend>
            <label className={backgroundMode === 'none' ? 'is-selected' : ''}><input type="radio" name="background" checked={backgroundMode === 'none'} onChange={() => { setBackgroundMode('none'); setResults(null) }} /><strong>No background</strong><small>Reports descriptive overlap metrics without P-values or FDR.</small></label>
            <label className={backgroundMode === 'union' ? 'is-selected' : ''}><input type="radio" name="background" checked={backgroundMode === 'union'} onChange={() => { setBackgroundMode('union'); setResults(null) }} /><strong>Input union</strong><small>Uses the observed union as a restricted statistical universe.</small></label>
            <label className={backgroundMode === 'upload' ? 'is-selected' : ''}><input type="radio" name="background" checked={backgroundMode === 'upload'} onChange={() => { setBackgroundMode('upload'); setResults(null) }} /><strong>Upload background</strong><small>Recommended. Include every site eligible for detection.</small></label>
          </fieldset>
          {mode === 'exact' && backgroundMode === 'upload' ? <label className="rnaenrich-background-upload"><Upload size={18} /><strong title={background?.name}>{background ? background.name : 'Choose background BED'}</strong><input type="file" accept=".bed,.txt,text/plain" onChange={(event) => { void uploadBackground(event.target.files?.[0]); event.target.value = '' }} /></label> : null}
          {mode === 'exact' && backgroundMode === 'union' ? <div className="rnaenrich-warning">Enrichment is evaluated only within the sites observed across the uploaded files. This restricted universe does not support transcriptome-wide inference.</div> : null}
          {mode === 'interval' ? <div className="rnaenrich-warning">Interval overlap reports descriptive effect sizes only. Valid permutation testing requires chromosome lengths and a callable-region mask, which cannot be inferred from the uploaded BED files.</div> : null}
          <div className="rnaenrich-action-row">
            <div className="rnaenrich-action-summary"><strong>{files.length} BED files</strong><span>{mode === 'exact' ? 'Exact match' : 'Interval overlap'}</span></div>
            <div className="rnaenrich-action-buttons">
              <button className="rnaenrich-button rnaenrich-button--quiet rnaenrich-demo-button" type="button" onClick={loadDemo}>Load demo</button>
              <button className="rnaenrich-button rnaenrich-button--primary" type="button" onClick={runAnalysis}>Run comparison</button>
            </div>
          </div>
          {message ? <div className="rnaenrich-message" role="status">{message}</div> : null}
        </div> : null}
      </section>

      {analysisTool === 'meta' ? <MetaPlot files={files} /> : null}

      {analysisTool === 'overlap' && results ? <section className="rnaenrich-results" aria-label="Comparison results">
        <article className="rnaenrich-panel">
          <div className="rnaenrich-panel-heading"><div><Grid3X3 size={18} /><h3>Pairwise similarity matrix</h3></div><div className="rnaenrich-panel-tools"><div className="rnaenrich-metric-tabs">
            <button className={metric === 'jaccard' ? 'is-active' : ''} type="button" onClick={() => setMetric('jaccard')}>Jaccard</button>
            <button className={metric === 'overlap' ? 'is-active' : ''} type="button" onClick={() => setMetric('overlap')}>Shared sites</button>
            <button className={metric === 'fdr' ? 'is-active' : ''} type="button" disabled={!results.some((pair) => pair.fdr !== null)} onClick={() => setMetric('fdr')}>−log10 FDR</button>
          </div><div className="rnaenrich-export-actions">
            <button type="button" onClick={() => exportFigure(buildMatrixSvg(files, results, metric), `rnaenrich_similarity_${metric}`, 'png')}><Download size={14} />PNG</button>
            <button type="button" onClick={() => exportFigure(buildMatrixSvg(files, results, metric), `rnaenrich_similarity_${metric}`, 'pdf')}><Download size={14} />PDF</button>
            <button type="button" onClick={exportMatrixCsv}><Download size={14} />CSV</button>
            <button type="button" onClick={exportMatrixTxt}><Download size={14} />TXT</button>
          </div></div></div>
          <div className="rnaenrich-export-surface"><Heatmap files={files} results={results} metric={metric} /></div>
        </article>

        {mode === 'exact' ? <article className="rnaenrich-panel">
          <div className="rnaenrich-panel-heading"><div><BarChart3 size={18} /><h3>Set intersections</h3></div><div className="rnaenrich-panel-tools"><span>{allCombinations.length > combinations.length ? `Top ${combinations.length} of ${allCombinations.length} combinations` : `${allCombinations.length} combination${allCombinations.length === 1 ? '' : 's'}`}</span><div className="rnaenrich-export-actions">
            <button type="button" onClick={() => exportFigure(buildUpSetSvg(files, combinations), 'rnaenrich_set_intersections', 'png')}><Download size={14} />PNG</button>
            <button type="button" onClick={() => exportFigure(buildUpSetSvg(files, combinations), 'rnaenrich_set_intersections', 'pdf')}><Download size={14} />PDF</button>
            <button type="button" onClick={exportIntersectionsCsv}><Download size={14} />CSV</button>
            <button type="button" onClick={exportIntersectionsTxt}><Download size={14} />TXT</button>
          </div></div></div>
          <div className="rnaenrich-export-surface"><UpSetPlot files={files} combinations={combinations} /></div>
        </article> : null}

        <article className="rnaenrich-panel">
          <div className="rnaenrich-panel-heading"><div><BarChart3 size={18} /><h3>Pairwise metrics</h3></div><div className="rnaenrich-export-actions">
            <button type="button" onClick={exportPairsCsv}><Download size={14} />CSV</button>
            <button type="button" onClick={exportPairsTxt}><Download size={14} />TXT</button>
          </div></div>
          <div className="rnaenrich-table-wrap"><table>
            <thead><tr><th>Pair</th><th>First sites</th><th>Second sites</th><th>Shared sites</th><th>First matched (%)</th><th>Second matched (%)</th><th>Jaccard</th><th>P-value</th><th>FDR</th><th>Shared BED</th></tr></thead>
            <tbody>{results.map((pair) => <tr key={`${pair.left}-${pair.right}`}>
              <td><strong title={pair.left}>{compactFilename(pair.left, files.findIndex((file) => file.name === pair.left))}</strong><span title={pair.right}>{compactFilename(pair.right, files.findIndex((file) => file.name === pair.right))}</span></td>
              <td>{pair.leftSize.toLocaleString()}</td><td>{pair.rightSize.toLocaleString()}</td><td>{pair.overlapCount.toLocaleString()}</td><td>{formatPercent(pair.leftRate)}</td><td>{formatPercent(pair.rightRate)}</td><td>{pair.jaccard.toFixed(4)}</td><td>{formatProbability(pair.pValue)}</td><td>{formatProbability(pair.fdr)}</td>
              <td>{mode === 'exact' ? <button type="button" onClick={() => {
                const left = files.find((file) => file.name === pair.left)
                const right = files.find((file) => file.name === pair.right)
                if (!left || !right) return
                const bed = exactSharedRecords(left, right).map((record) => `${record.chrom}\t${record.start}\t${record.end}`).join('\n')
                downloadText(`${safeFilename(pair.left)}_${safeFilename(pair.right)}_shared.bed`, bed ? `${bed}\n` : '')
              }}><Download size={14} /> BED</button> : '—'}</td>
            </tr>)}</tbody>
          </table></div>
        </article>
      </section> : null}
    </main>

    <footer className="rnaenrich-footer">RNAmeta | Rice Research Institute, Guangdong Academy of Agricultural Sciences | Guangzhou, Guangdong, China | 2026</footer>
  </div>
}
