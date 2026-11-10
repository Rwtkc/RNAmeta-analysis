import { Download, FileUp, ExternalLink } from 'lucide-react'
import { jsPDF } from 'jspdf'
import { svg2pdf } from 'svg2pdf.js'
import { useEffect, useMemo, useRef, useState } from 'react'
import './RNAPlot.css'
import { AnnotationTracks, StructureViewer } from './rnaAnnotationMigration'
import { buildAnnotationLegendMarkup, type StructureHighlight } from './rnaAnnotationUtils'

const appBase = import.meta.env.BASE_URL

const DEMO_SEQUENCE = 'AUGGCUACGUACGGAUCUAGCUAUGCAUGGAUACGGCUAGUACCGUAUGCAAGCUUGACUAGGCUAACGUAGCUAUGGACUACGUAGCUAGCAUCGUACGGAUACUAGCUACGUA'
const DEMO_ANNOTATIONS = `TR1 1 2 m6A
TR1 2 3 psi
TR1 4 5 m5C
TR1 6 7 m1A
TR1 9 10 psi
TR1 10 12 m1A
TR1 11 14 m6A
TR1 12 15 m6A peak
TR1 13 14 m5C
TR1 14 16 psi
TR1 17 19 m1A
TR1 18 20 m5C
TR1 22 23 m6A
TR1 26 27 m6A
TR1 27 29 m5C
TR1 35 36 psi
TR1 39 40 m1A
TR1 42 45 psi
TR1 44 48 m6A peak
TR1 46 47 m5C
TR1 47 50 psi
TR1 48 52 m6A
TR1 49 51 m6A
TR1 50 53 m5C
TR1 52 54 m1A
TR1 53 55 psi
TR1 57 58 m5C
TR1 60 62 m1A
TR1 65 66 m6A
TR1 68 69 m5C
TR1 72 74 m1A
TR1 73 76 m5C
TR1 74 75 psi
TR1 75 78 m6A
TR1 76 79 m6A
TR1 77 78 m1A
TR1 79 82 m6A peak
TR1 80 84 psi
TR1 88 89 m5C
TR1 90 91 psi
TR1 92 95 m1A
TR1 97 98 m6A
TR1 99 100 m1A
TR1 101 104 m5C
TR1 103 106 m6A peak
TR1 104 105 m5C
TR1 105 107 psi
TR1 108 110 m1A
TR1 111 112 m6A`
const DEMO_REGIONS = `TR1 0 16 5' UTR
TR1 20 91 CDS
TR1 96 115 3' UTR`

type Mode = 'sequence' | 'annotations'
type AnnotationRow = { transcript: string; start: number; end: number; label: string }
type RegionRow = AnnotationRow
type Result = {
  sequence: string
  annotations: AnnotationRow[]
  regions: RegionRow[]
  plot_data_uri?: string
  structure_payload?: { svg_data_uri?: string; svg_markup?: string; method?: string; sequence?: string; highlights?: StructureHighlight[] }
  summary?: Record<string, unknown>
}

function normalizeSequence(value: string) {
  return value.split(/\r?\n/).filter((line) => !line.trim().startsWith('>')).join('').replace(/\s+/g, '').toUpperCase().replaceAll('T', 'U')
}

function decodeSvgDataUri(uri: string) {
  const separator = uri.indexOf(',')
  if (separator < 0) return null
  const metadata = uri.slice(0, separator)
  const payload = uri.slice(separator + 1)
  if (!metadata.includes('image/svg+xml')) return null
  if (!metadata.includes(';base64')) return decodeURIComponent(payload)
  const binary = atob(payload)
  const bytes = Uint8Array.from(binary, (character) => character.charCodeAt(0))
  return new TextDecoder().decode(bytes)
}

function encodeSvgDataUri(svg: string) {
  const bytes = new TextEncoder().encode(svg)
  let binary = ''
  bytes.forEach((byte) => { binary += String.fromCharCode(byte) })
  return `data:image/svg+xml;base64,${btoa(binary)}`
}

const PLOT_ANNOTATION_COLORS: Record<string, string> = {
  m6a: '#7fc97f',
  psi: '#e6a0a0',
  m5c: '#b39ddb',
  m1a: '#80cbc4',
  'm6a peak': '#81c784'
}

const PLOT_REGION_COLORS: Record<string, string> = {
  "5' utr": '#4a90b5',
  cds: '#4ca34c',
  "3' utr": '#c75b60'
}

function appendPlotLegend(uri: string, annotations: AnnotationRow[], regions: RegionRow[]) {
  const source = decodeSvgDataUri(uri)
  if (!source) return uri
  const documentNode = new DOMParser().parseFromString(source, 'image/svg+xml')
  const svg = documentNode.documentElement
  if (!svg || svg.nodeName.toLowerCase() !== 'svg') return uri
  svg.querySelectorAll('rect').forEach((rect) => {
    if (rect.getAttribute('class')) return
    const style = (rect.getAttribute('style') || '').toLowerCase().replace(/\s+/g, '')
    const fill = (rect.getAttribute('fill') || '').toLowerCase()
    const isWhite = style.includes('fill:#ffffff') || style.includes('fill:white') || fill === '#ffffff' || fill === 'white'
    const hasStroke = style.includes('stroke:#ffffff') || style.includes('stroke:white') || (rect.getAttribute('stroke') || '').toLowerCase() === '#ffffff'
    const rectX = (rect.getAttribute('x') || '').trim()
    const rectY = (rect.getAttribute('y') || '').trim()
    const rectWidth = (rect.getAttribute('width') || '').trim()
    const rectHeight = (rect.getAttribute('height') || '').trim()
    const isCanvasRect = (rectX === '' || rectX === '0' || rectX === '0.00') && (rectY === '' || rectY === '0' || rectY === '0.00') && (rectWidth === '100%' || rectHeight === '100%' || (Number.parseFloat(rectWidth) >= 800 && Number.parseFloat(rectHeight) >= 300))
    if (isWhite && (isCanvasRect || !hasStroke || style.includes('stroke:none'))) rect.remove()
  })
  const labels = [...new Set(annotations.map((row) => row.label.trim()))]
  const annotationItems = labels.map((label, index) => ({ label, color: PLOT_ANNOTATION_COLORS[label.toLowerCase()] || ['#7fc97f', '#e6a0a0', '#b39ddb', '#80cbc4', '#81c784'][index % 5] }))
  const regionItems = [...new Set(regions.map((row) => row.label.trim()))].map((label, index) => ({ label, color: PLOT_REGION_COLORS[label.toLowerCase()] || ['#4a90b5', '#4ca34c', '#c75b60'][index % 3] }))
  const viewBox = (svg.getAttribute('viewBox') || '').trim().split(/[\s,]+/).map(Number)
  const [viewX, viewY, viewWidth, viewHeight] = viewBox.length === 4 && viewBox.every(Number.isFinite)
    ? viewBox
    : [0, 0, Number.parseFloat(svg.getAttribute('width') || '864'), Number.parseFloat(svg.getAttribute('height') || '331')]
  const legendItems = [...annotationItems, ...regionItems]
  const legendWidth = Math.min(700, Math.max(240, viewWidth - 28))
  const legendColumns = Math.max(1, Math.floor(legendWidth / 72))
  const legendRows = Math.ceil(legendItems.length / legendColumns)
  const legendHeight = Math.max(32, legendRows * 15 + 12)
  const legend = buildAnnotationLegendMarkup(legendItems, viewX + (viewWidth - legendWidth) / 2, viewY + viewHeight + 15, legendWidth, 'center')
  const viewBoxHeight = viewHeight + legendHeight
  svg.setAttribute('viewBox', `${viewX} ${viewY} ${viewWidth} ${viewBoxHeight}`)
  const svgHeight = svg.getAttribute('height') || ''
  const heightUnit = svgHeight.match(/[a-z%]+$/i)?.[0] || ''
  svg.setAttribute('height', `${viewBoxHeight}${heightUnit}`)
  if (legend) {
    const legendDocument = new DOMParser().parseFromString(`<svg xmlns="http://www.w3.org/2000/svg">${legend}</svg>`, 'image/svg+xml')
    const legendNode = legendDocument.documentElement.firstElementChild
    if (legendNode) svg.appendChild(documentNode.importNode(legendNode, true))
  }
  return encodeSvgDataUri(new XMLSerializer().serializeToString(svg))
}

function parseRows(value: string): AnnotationRow[] {
  return value.split(/\r?\n/).map((line) => line.trim()).filter((line) => line && !line.startsWith('#')).flatMap((line) => {
    const fields = line.split(/\s+/)
    const start = Number.parseInt(fields[1] || '', 10)
    const end = Number.parseInt(fields[2] || '', 10)
    const label = fields.slice(3).join(' ').trim()
    return fields[0] && Number.isFinite(start) && Number.isFinite(end) && label ? [{ transcript: fields[0], start, end, label }] : []
  })
}

/* RNAmeta structure rendering is implemented in rnaAnnotationMigration.tsx. */

function triggerDownload(blob: Blob, filename: string) {
  const url = URL.createObjectURL(blob)
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = filename
  document.body.appendChild(anchor)
  anchor.click()
  anchor.remove()
  window.setTimeout(() => URL.revokeObjectURL(url), 0)
}

const EXPORT_STYLE_PROPERTIES = ['fill', 'fill-opacity', 'stroke', 'stroke-opacity', 'stroke-width', 'stroke-linecap', 'stroke-linejoin', 'stroke-dasharray', 'stroke-dashoffset', 'opacity', 'font-family', 'font-size', 'font-weight', 'font-style', 'text-anchor', 'dominant-baseline', 'alignment-baseline', 'paint-order', 'display', 'visibility', 'shape-rendering']

function cloneSvgForExport(node: SVGSVGElement) {
  const clone = node.cloneNode(true) as SVGSVGElement
  const sourceElements = [node, ...Array.from(node.querySelectorAll('*'))]
  const targetElements = [clone, ...Array.from(clone.querySelectorAll('*'))]
  sourceElements.forEach((source, index) => {
    const target = targetElements[index]
    if (!(source instanceof SVGElement) || !(target instanceof SVGElement)) return
    const computed = getComputedStyle(source)
    EXPORT_STYLE_PROPERTIES.forEach((property) => {
      const value = computed.getPropertyValue(property)
      if (value) target.style.setProperty(property, value)
    })
  })
  return clone
}

function prepareSvgForExport(svg: SVGSVGElement) {
  svg.querySelectorAll('.rna-annotation-tracks__band').forEach((element) => element.setAttribute('fill', '#f8f8f3'))
  svg.querySelectorAll('.rna-annotation-tracks__grid').forEach((element) => { element.setAttribute('stroke', '#e4e6df'); element.setAttribute('stroke-width', '0.75') })
  svg.querySelectorAll('.rna-annotation-tracks__baseline').forEach((element) => { element.setAttribute('stroke', '#d2d1cb'); element.setAttribute('stroke-width', '0.9') })
  svg.querySelectorAll('.rna-annotation-tracks__axis, .rna-annotation-tracks__tick').forEach((element) => { element.setAttribute('stroke', '#667067'); element.setAttribute('stroke-width', '1.1') })
  svg.querySelectorAll('.rna-annotation-tracks__interval').forEach((element) => {
    element.setAttribute('fill', '#243229')
    element.setAttribute('stroke', 'none')
    element.setAttribute('stroke-width', '0')
    element.setAttribute('font-family', 'Arial, Helvetica, sans-serif')
    element.setAttribute('font-size', '9px')
    element.setAttribute('font-weight', '400')
    ;(element as SVGElement).style.setProperty('stroke', 'none')
    ;(element as SVGElement).style.setProperty('stroke-width', '0')
  })
  svg.querySelectorAll('.rnameta-marked-base-label, .numbering-label.sequential').forEach((element) => {
    element.setAttribute('fill', '#000000')
    element.setAttribute('stroke', 'none')
    element.setAttribute('stroke-width', '0')
    element.setAttribute('font-weight', '400')
    ;(element as SVGElement).style.setProperty('fill', '#000000')
    ;(element as SVGElement).style.setProperty('stroke', 'none')
    ;(element as SVGElement).style.setProperty('stroke-width', '0')
    ;(element as SVGElement).style.setProperty('font-weight', '400')
  })
  svg.querySelectorAll('.rnameta-position-labels text, .rnameta-terminal-labels text, text').forEach((element) => {
    if (!element.getAttribute('font-weight')) element.setAttribute('font-weight', '400')
  })
  return svg
}

function svgToPngBlob(node: SVGSVGElement): Promise<Blob> {
  const { width, height } = getSvgDimensions(node)
  const clone = prepareSvgForExport(cloneSvgForExport(node))
  clone.setAttribute('width', String(width))
  clone.setAttribute('height', String(height))
  clone.style.width = `${width}px`
  clone.style.height = `${height}px`
  const source = new XMLSerializer().serializeToString(clone)
  const url = URL.createObjectURL(new Blob([source], { type: 'image/svg+xml;charset=utf-8' }))
  const image = new Image()
  return new Promise((resolve, reject) => {
    image.onerror = () => { URL.revokeObjectURL(url); reject(new Error('Unable to render SVG as PNG.')) }
    image.onload = () => {
      const canvas = document.createElement('canvas')
      const scale = 3
      canvas.width = Math.ceil(width * scale)
      canvas.height = Math.ceil(height * scale)
      const context = canvas.getContext('2d')
      if (!context) { URL.revokeObjectURL(url); reject(new Error('Unable to create PNG canvas.')); return }
      context.fillStyle = '#fffdf8'
      context.fillRect(0, 0, canvas.width, canvas.height)
      context.scale(scale, scale)
      context.drawImage(image, 0, 0, width, height)
      canvas.toBlob((blob) => {
        URL.revokeObjectURL(url)
        if (!blob) { reject(new Error('Unable to encode PNG.')); return }
        resolve(blob)
      }, 'image/png')
    }
    image.src = url
  })
}

function getSvgDimensions(node: SVGSVGElement) {
  const viewBox = node.viewBox.baseVal
  const width = viewBox.width || Number.parseFloat(node.getAttribute('width') || '')
  const height = viewBox.height || Number.parseFloat(node.getAttribute('height') || '')
  if (!(width > 0 && height > 0)) throw new Error('Unable to determine SVG dimensions.')
  return { width, height }
}

function isRenderedSvg(element: Element | null): element is SVGSVGElement {
  return Boolean(element && element.localName === 'svg' && element.namespaceURI === 'http://www.w3.org/2000/svg')
}

async function svgToPdfBlob(node: SVGSVGElement) {
  const { width, height } = getSvgDimensions(node)
  const clone = prepareSvgForExport(cloneSvgForExport(node))
  clone.setAttribute('width', String(width))
  clone.setAttribute('height', String(height))
  clone.style.width = `${width}px`
  clone.style.height = `${height}px`
  clone.setAttribute('font-family', 'Helvetica, Arial, sans-serif')
  clone.querySelectorAll('text, tspan').forEach((text) => {
    text.setAttribute('font-family', 'Helvetica, Arial, sans-serif')
    const textElement = text as SVGElement
    textElement.style.fontFamily = 'Helvetica, Arial, sans-serif'
  })
  const pdf = new jsPDF({ orientation: width >= height ? 'landscape' : 'portrait', unit: 'px', format: [width, height], compress: true })
  await svg2pdf(clone, pdf, { x: 0, y: 0, width, height })
  return pdf.output('blob')
}

async function imageSourceToSvg(source: string) {
  const response = await fetch(source)
  if (!response.ok) throw new Error('Unable to read the annotation graphic.')
  const text = await response.text()
  const svg = new DOMParser().parseFromString(text, 'image/svg+xml').documentElement
  return svg.nodeName.toLowerCase() === 'svg' ? svg as unknown as SVGSVGElement : null
}

function parseInputResult(mode: Mode, sequence: string, annotations: string, regions: string): Result {
  const normalized = mode === 'sequence' ? normalizeSequence(sequence) : ''
  const marks = parseRows(annotations)
  const regionRows = parseRows(regions)
  const length = normalized.length || Math.max(1, ...marks.map((row) => row.end), ...regionRows.map((row) => row.end))
  return { sequence: normalized || 'N'.repeat(length), annotations: marks, regions: regionRows }
}

export default function RNAPlot() {
  const [mode, setMode] = useState<Mode>('sequence')
  const [sequence, setSequence] = useState('')
  const [sequenceAnnotations, setSequenceAnnotations] = useState('')
  const [annotations, setAnnotations] = useState('')
  const [regions, setRegions] = useState('')
  const [results, setResults] = useState<Record<Mode, Result | null>>({ sequence: null, annotations: null })
  const [busy, setBusy] = useState(false)
  const [notice, setNotice] = useState('')
  const outputRef = useRef<HTMLElement>(null)
  const fieldsRef = useRef<HTMLDivElement>(null)
  const activeTextareaRef = useRef<HTMLTextAreaElement | null>(null)
  const activeAnnotations = mode === 'sequence' ? sequenceAnnotations : annotations
  const result = results[mode]
  const metrics = useMemo(() => {
    const rows = parseRows(activeAnnotations)
    const regionRows = mode === 'annotations' ? parseRows(regions) : []
    return { transcript: rows[0]?.transcript || regionRows[0]?.transcript || '', length: mode === 'sequence' ? normalizeSequence(sequence).length : Math.max(0, ...rows.map((row) => row.end), ...regionRows.map((row) => row.end)), marks: rows.length, types: new Set(rows.map((row) => row.label.toLowerCase())).size }
  }, [activeAnnotations, mode, regions, sequence])

  useEffect(() => {
    const container = fieldsRef.current
    if (!container) return
    const syncHeights = () => {
      const textareas = Array.from(container.querySelectorAll('textarea'))
      if (!textareas.length) return
      const active = activeTextareaRef.current && textareas.includes(activeTextareaRef.current) ? activeTextareaRef.current : textareas[0]
      const height = active.getBoundingClientRect().height
      textareas.forEach((textarea) => {
        if (textarea !== active) textarea.style.height = `${height}px`
      })
    }
    const observer = new ResizeObserver(syncHeights)
    container.querySelectorAll('textarea').forEach((textarea) => {
      textarea.addEventListener('pointerdown', () => { activeTextareaRef.current = textarea })
      observer.observe(textarea)
    })
    syncHeights()
    return () => observer.disconnect()
  }, [mode])

  const renderPreview = async () => {
    setBusy(true)
    setNotice('')
    const payload = { input_mode: mode, transcript_id: metrics.transcript || 'TR1', sequence: mode === 'sequence' ? sequence : '', annotations: activeAnnotations, regions: mode === 'annotations' ? regions : '' }
    try {
      const controller = new AbortController()
      const timer = window.setTimeout(() => controller.abort(), 120000)
      const response = await fetch(`${appBase}RNAplot/api/render`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload), signal: controller.signal })
      window.clearTimeout(timer)
      const remote = await response.json() as Partial<Result> & { error?: string }
      if (!response.ok) throw new Error(remote.error || 'RNAplot calculation failed.')
      if (!remote.summary) throw new Error('RNAplot returned an incomplete calculation result.')
      const parsedInput = parseInputResult(mode, sequence, activeAnnotations, regions)
      const normalizedRemote = remote.plot_data_uri ? { ...remote, plot_data_uri: appendPlotLegend(remote.plot_data_uri, parsedInput.annotations, parsedInput.regions) } : remote
      setResults((previous) => ({ ...previous, [mode]: { ...parsedInput, ...normalizedRemote } }))
    } catch (error) {
      setResults((previous) => ({ ...previous, [mode]: null }))
      setNotice(error instanceof Error ? error.message : 'RNAplot calculation failed.')
    } finally {
      setBusy(false)
    }
  }

  const loadDemo = () => {
    if (mode === 'sequence') { setSequence(DEMO_SEQUENCE); setSequenceAnnotations(DEMO_ANNOTATIONS) }
    else { setAnnotations(DEMO_ANNOTATIONS); setRegions(DEMO_REGIONS) }
    setResults((previous) => ({ ...previous, [mode]: null }))
    setNotice('')
  }

  const changeMode = (nextMode: Mode) => {
    if (nextMode === mode) return
    setMode(nextMode)
    setNotice('')
  }

  const importText = (event: React.ChangeEvent<HTMLInputElement>, setter: (value: string) => void) => {
    const file = event.target.files?.[0]
    event.target.value = ''
    if (file) file.text().then(setter).catch(() => setNotice('无法读取该文件。'))
  }

  const exportOutput = async (format: 'png' | 'pdf') => {
    const transcript = (metrics.transcript || 'TR1').replace(/[^A-Za-z0-9._-]+/g, '_')
    try {
      setNotice('')
      const output = outputRef.current
      if (!output) return
      const toBlob = (node: SVGSVGElement) => format === 'png' ? svgToPngBlob(node) : svgToPdfBlob(node)
      if (mode === 'sequence') {
        // The tracks are rendered as SVG; validate the namespace explicitly so
        // the export also works when the node came from dynamic markup.
        const tracks = output.querySelector('.rna-annotation-tracks svg')
        if (!isRenderedSvg(tracks)) throw new Error('The annotation tracks must finish rendering before export.')
        triggerDownload(await toBlob(tracks), `rnaplot-${transcript}-annotation-tracks.${format}`)
        return
      }
      const image = output.querySelector('.rna-annotation-static-plot img')
      if (!(image instanceof HTMLImageElement)) throw new Error('The annotation map must finish rendering before export.')
      const source = image.currentSrc || image.src
      const svg = await imageSourceToSvg(source)
      if (svg) {
        triggerDownload(await toBlob(svg), `rnaplot-${transcript}-annotation-map.${format}`)
        return
      }
      if (format === 'png') {
        const response = await fetch(source)
        if (!response.ok) throw new Error('Unable to read the annotation map.')
        triggerDownload(await response.blob(), `rnaplot-${transcript}-annotation-map.png`)
        return
      }
      throw new Error('The annotation map is not an SVG and cannot be exported as a vector PDF.')
    } catch (error) {
      setNotice(error instanceof Error ? error.message : 'Export failed.')
    }
  }

  return <div className="rnaplot-page">
    <header className="rnaplot-header">
      <a className="rnaplot-brand" href={appBase}>RNA<span>meta</span></a>
      <nav className="rnaplot-nav" aria-label="Main navigation">
        <a href={appBase}>Home</a>
        <a href={`${appBase}#news`}>News</a>
        <a href="https://github.com/Rwtkc/RNAmeta_R" target="_blank" rel="noreferrer">GitHub <ExternalLink size={13} aria-hidden="true" /></a>
      </nav>
    </header>
    <main className="rnaplot-main">
      <section className="rnaplot-heading" aria-label="About RNAplot">
        <div className="rnaplot-intro-copy"><strong>RNAplot helps you turn transcript-relative RNA modification coordinates into clear visual maps.</strong><p>Use <em>RNA sequence</em> mode when you have an RNA sequence and a list of annotated positions. Use <em>Annotations</em> mode when you only have transcript-relative intervals, with optional transcript regions such as 5' UTR, CDS, or 3' UTR. Enter one row per feature in the format <code>TR1 start end label</code>, then choose <strong>Render preview</strong> to generate the annotation tracks and structure context.</p></div>
      </section>

      <section className="rnaplot-panel rnaplot-panel--input" aria-labelledby="rnaplot-input-title">
        <div className="rnaplot-section-heading">
          <h2 id="rnaplot-input-title">Prepare your annotation input</h2>
          <span className="rnaplot-section-label">Input</span>
        </div>
        <div className="rnaplot-input-panel">
          <span className="rnaplot-mode-label">Select input type</span>
          <div className="rnaplot-mode-switch" role="tablist"><button type="button" role="tab" aria-selected={mode === 'sequence'} className={mode === 'sequence' ? 'is-active' : ''} onClick={() => changeMode('sequence')}>RNA sequence</button><button type="button" role="tab" aria-selected={mode === 'annotations'} className={mode === 'annotations' ? 'is-active' : ''} onClick={() => changeMode('annotations')}>Annotations</button></div>
          <div ref={fieldsRef} className={`rnaplot-fields rnaplot-fields--${mode}`}>
            {mode === 'sequence' ? <>
              <label className="rnaplot-field"><span>RNA sequence</span><textarea value={sequence} onChange={(event) => setSequence(event.target.value)} spellCheck="false" placeholder="Paste FASTA or an RNA sequence" /></label>
              <label className="rnaplot-field"><span>Annotations</span><textarea value={sequenceAnnotations} onChange={(event) => setSequenceAnnotations(event.target.value)} spellCheck="false" placeholder="TR1 12 15 m6A" /></label>
            </> : <label className="rnaplot-field"><span>Annotations</span><textarea value={annotations} onChange={(event) => setAnnotations(event.target.value)} spellCheck="false" placeholder="TR1 12 15 m6A" /></label>}
            {mode === 'annotations' ? <label className="rnaplot-field"><span>Transcript regions</span><textarea className="rnaplot-region-input" value={regions} onChange={(event) => setRegions(event.target.value)} spellCheck="false" placeholder={"TR1 0 16 5' UTR"} /></label> : null}
          </div>
          <input className="rnaplot-hidden-file" type="file" accept=".txt,text/plain" id="rnaplot-upload" onChange={(event) => importText(event, mode === 'sequence' ? setSequenceAnnotations : setAnnotations)} />
          <div className="rnaplot-actions"><label className="rnaplot-button rnaplot-button--soft" htmlFor="rnaplot-upload"><FileUp size={16} />Upload annotations</label><button type="button" className="rnaplot-button rnaplot-button--primary" onClick={renderPreview} disabled={busy}><span>{busy ? 'Rendering...' : 'Render preview'}</span></button><button type="button" className="rnaplot-button rnaplot-button--quiet" onClick={loadDemo} disabled={busy}>Load demo</button></div>
        </div>
      </section>

      {result ? <section ref={outputRef} className="rnaplot-panel rnaplot-panel--output" aria-label="Annotation output">
        <div className="rnaplot-output-topline"><div className="rnaplot-output-tools"><button className="rnaplot-export-button" type="button" title="导出 PNG" aria-label="导出 PNG" onClick={() => void exportOutput('png')}><Download size={15} />PNG</button><button className="rnaplot-export-button" type="button" title="导出 PDF" aria-label="导出 PDF" onClick={() => void exportOutput('pdf')}><Download size={15} />PDF</button></div></div>
        {notice ? <div className="rnaplot-notice">{notice}</div> : null}
        {mode === 'sequence' ? <div className="rna-annotation-structure-preview"><div className="rna-annotation-figure__canvas">{result.structure_payload?.svg_data_uri ? <StructureViewer src={result.structure_payload.svg_data_uri} markup={result.structure_payload.svg_markup} highlights={result.structure_payload.highlights || []} /> : null}</div><div className="rna-annotation-figure--tracks"><AnnotationTracks highlights={result.structure_payload?.highlights || []} sequenceLength={String(result.structure_payload?.sequence || result.sequence).length || result.sequence.length} /></div></div> : <div className="rna-annotation-static-plot"><div className="rna-annotation-static-plot__canvas">{result.plot_data_uri ? <img className="rnaplot-remote-plot" src={result.plot_data_uri} alt="RNA annotation map" /> : null}</div></div>}
      </section> : null}
    </main>
    <footer className="rnaplot-footer">RNAmeta | Rice Research Institute, Guangdong Academy of Agricultural Sciences | Guangzhou, Guangdong, China | 2026</footer>
  </div>
}
