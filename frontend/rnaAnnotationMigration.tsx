import { useEffect, useMemo, useRef, useState } from 'react'

type StructureHighlight = { start?: number; end?: number; label?: string; category?: string }

const ANNOTATION_COLORS = ['#525C81', '#A6DCDD', '#AD5490', '#A45C5F', '#E9C064', '#678562']
const ANNOTATION_COLOR_BY_LABEL: Record<string, string> = {
  m6a: '#525C81', m5c: '#A6DCDD', m1a: '#AD5490', psi: '#A45C5F', pseudouridine: '#A45C5F', 'm6a peak': '#E9C064'
}

function parseStructureCoordinates(svgRoot: SVGSVGElement | Element) {
  const titledCoordinates = Array.from(svgRoot.querySelectorAll('g')).flatMap((group) => {
    const title = group.querySelector(':scope > title')?.textContent || ''
    const match = title.match(/^(\d+) \(position\.label/)
    const text = group.querySelector(':scope > text')
    if (!match || !text) return []
    const nucleotide = text.textContent?.trim().toUpperCase() || ''
    if (!/^[ACGUTN]$/.test(nucleotide)) return []
    const x = Number.parseFloat(text.getAttribute('x') || '')
    const y = Number.parseFloat(text.getAttribute('y') || '')
    return Number.isFinite(x) && Number.isFinite(y) ? [{ position: Number(match[1]), x, y, nucleotide }] : []
  })
  if (titledCoordinates.length) return titledCoordinates

  // Current R2DT/RNApuzzler output uses one transformed #seq group with
  // sequential text.nucleotide nodes instead of titled position groups.
  return Array.from(svgRoot.querySelectorAll('text.nucleotide')).flatMap((text, index) => {
    const nucleotide = text.textContent?.trim().toUpperCase() || ''
    const x = Number.parseFloat(text.getAttribute('x') || '')
    const y = Number.parseFloat(text.getAttribute('y') || '')
    return /^[ACGUTN]$/.test(nucleotide) && Number.isFinite(x) && Number.isFinite(y)
      ? [{ position: index + 1, x, y, nucleotide }]
      : []
  })
}

function darkerColor(color: string) {
  const match = String(color || '').match(/^#([0-9a-f]{6})$/i)
  if (!match) return '#50504c'
  const values = [0, 2, 4].map((offset) => Number.parseInt(match[1].slice(offset, offset + 2), 16))
  return `#${values.map((value) => Math.round(value * .64).toString(16).padStart(2, '0')).join('')}`
}

function hslToHex(hue: number, saturation: number, lightness: number) {
  const chroma = (1 - Math.abs(2 * lightness - 1)) * saturation
  const segment = hue / 60
  const secondary = chroma * (1 - Math.abs((segment % 2) - 1))
  const [red, green, blue] = segment < 1 ? [chroma, secondary, 0] : segment < 2 ? [secondary, chroma, 0] : segment < 3 ? [0, chroma, secondary] : segment < 4 ? [0, secondary, chroma] : segment < 5 ? [secondary, 0, chroma] : [chroma, 0, secondary]
  const offset = lightness - chroma / 2
  return `#${[red, green, blue].map((value) => Math.round((value + offset) * 255).toString(16).padStart(2, '0')).join('')}`
}

function labelHash(value: string) { return Array.from(value).reduce((hash, character) => ((hash * 31) + character.charCodeAt(0)) >>> 0, 7) }

function colorForLabel(label: string, usedColors: Set<string>) {
  if (ANNOTATION_COLOR_BY_LABEL[label]) return ANNOTATION_COLOR_BY_LABEL[label]
  const seed = labelHash(label)
  for (let attempt = 0; attempt < 36; attempt += 1) {
    const hue = (seed + attempt * 137.508) % 360
    const color = hslToHex(hue, 0.47 + ((seed >>> 4) % 10) / 100, 0.62 + ((seed >>> 8) % 8) / 100)
    if (!usedColors.has(color)) return color
  }
  return ANNOTATION_COLORS[seed % ANNOTATION_COLORS.length]
}

function escapeXml(value: string) {
  return value.replace(/[&<>"']/g, (character) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&apos;' }[character] || character))
}

function buildAnnotationLegendMarkup(items: Array<{ label: string; color: string }>, x: number, y: number, maxWidth = 180, align: 'start' | 'center' = 'start') {
  if (!items.length) return ''
  const itemWidth = Math.max(48, ...items.map((item) => item.label.length * 4.8 + 22))
  const rowHeight = 15
  let row = 0
  let column = 0
  const parts = items.map((item) => {
    if (column > 0 && (column + 1) * itemWidth > maxWidth) { row += 1; column = 0 }
    const columns = Math.max(1, Math.floor(maxWidth / itemWidth))
    const rowStart = row * columns
    const rowCount = Math.min(columns, items.length - rowStart)
    const rowOffset = align === 'center' ? Math.max(0, (maxWidth - rowCount * itemWidth) / 2) : 0
    const itemX = x + rowOffset + column * itemWidth
    const itemY = y + row * rowHeight
    column += 1
    return `<g class="rnameta-legend-item"><circle cx="${itemX.toFixed(2)}" cy="${itemY.toFixed(2)}" r="3.7" fill="${item.color}" fill-opacity="1" stroke="${darkerColor(item.color)}" stroke-width="0.75" style="fill:${item.color} !important;fill-opacity:1 !important;stroke:${darkerColor(item.color)} !important;stroke-width:0.75 !important"/><text x="${(itemX + 7).toFixed(2)}" y="${(itemY + 1.9).toFixed(2)}" fill="#243229" font-family="Arial, Helvetica, sans-serif" font-size="6.2" font-weight="400" style="fill:#243229 !important;font-weight:400">${escapeXml(item.label)}</text></g>`
  }).join('')
  return `<g class="rnameta-legend" pointer-events="none">${parts}</g>`
}

function buildStructureMarkup(markup: string, highlights: StructureHighlight[] = []) {
  if (!markup) return ''
  const documentNode = new DOMParser().parseFromString(markup, 'image/svg+xml')
  const svg = documentNode.documentElement
  if (!svg || svg.nodeName.toLowerCase() !== 'svg') return markup
  // R2DT emits a full-canvas white rect for standalone image export. RNAmeta
  // displays the structure on its own transparent viewer surface, so remove
  // only that canvas background and keep any structural rectangles intact.
  const rootWidth = Number.parseFloat(svg.getAttribute('width') || '')
  const rootHeight = Number.parseFloat(svg.getAttribute('height') || '')
  svg.querySelectorAll('rect').forEach((rect) => {
    const style = (rect.getAttribute('style') || '').toLowerCase().replace(/\s+/g, '')
    const isFullCanvas = rect.getAttribute('x') === '0' && rect.getAttribute('y') === '0' && Number.parseFloat(rect.getAttribute('width') || '') === rootWidth && Number.parseFloat(rect.getAttribute('height') || '') === rootHeight
    const isWhiteCanvas = style.includes('fill:white') && style.includes('stroke:white') && (!rect.getAttribute('class'))
    if (isFullCanvas || isWhiteCanvas) rect.remove()
  })
  const coordinates = parseStructureCoordinates(svg)
  if (!coordinates.length) return markup
  const hasTitledCoordinates = Array.from(svg.querySelectorAll('g')).some((group) => /^\d+ \(position\.label/.test(group.querySelector(':scope > title')?.textContent || ''))
  const nucleotideTexts = Array.from(svg.querySelectorAll('text.nucleotide'))
  const nucleotideParent = nucleotideTexts[0]?.parentElement || null
  const renderParent = hasTitledCoordinates ? svg : (nucleotideParent?.parentElement || svg)
  // RNAmeta renders its annotation overlays in the untransformed SVG
  // coordinate system. New R2DT output wraps the drawing in a scale/translate
  // group, so remove that wrapper before adding the same overlays; this keeps
  // line geometry, circle radii, and text sizes in one coordinate system.
  if (!hasTitledCoordinates && renderParent !== svg) renderParent.removeAttribute('transform')
  const renderCoordinates = coordinates
  const labels = [...new Set(highlights.map((annotation) => String(annotation.label || annotation.category || 'annotation').trim().toLowerCase()))].sort()
  const usedColors = new Set<string>()
  const colorByLabel = new Map(labels.map((label) => { const color = colorForLabel(label, usedColors); usedColors.add(color); return [label, color] as const }))
  const legendItems = labels.map((label) => ({ label, color: colorByLabel.get(label) || ANNOTATION_COLORS[0] }))
  const markedPositions = new Map<number, string>()
  highlights.forEach((annotation) => {
    const key = String(annotation.label || annotation.category || 'annotation').trim().toLowerCase()
    const color = colorByLabel.get(key) || ANNOTATION_COLORS[0]
    for (let position = Number(annotation.start); position <= Number(annotation.end); position += 1) markedPositions.set(position, color)
  })
  const baseNodes = renderCoordinates.map((point) => { const fill = markedPositions.get(point.position) || '#d8d8d3'; const stroke = markedPositions.has(point.position) ? darkerColor(fill) : '#767671'; return `<circle class="rnameta-base-node" cx="${point.x.toFixed(3)}" cy="${point.y.toFixed(3)}" r="6.2" fill="${fill}" fill-opacity="1" stroke="${stroke}" stroke-width="1.15"/>` }).join('')
  const nativeNumberLabels = new Set(Array.from(svg.querySelectorAll('text.numbering-label')).map((text) => Number.parseInt(text.textContent?.trim() || '', 10)).filter(Number.isFinite))
  const center = renderCoordinates.reduce((total, point) => ({ x: total.x + point.x, y: total.y + point.y }), { x: 0, y: 0 }); center.x /= renderCoordinates.length; center.y /= renderCoordinates.length
  const markedLetters = renderCoordinates.filter((point) => markedPositions.has(point.position)).map((point) => `<text class="rnameta-marked-base-label" x="${point.x.toFixed(3)}" y="${point.y.toFixed(3)}" fill="#000000" stroke="#ffffff" stroke-width="1.1" paint-order="stroke fill" font-size="4.4" font-family="Arial, Helvetica, sans-serif" font-weight="400" text-anchor="middle" dominant-baseline="middle" alignment-baseline="middle" style="font-weight:400;fill:#000000;stroke:#ffffff;stroke-width:1.1;paint-order:stroke fill;text-anchor:middle;dominant-baseline:middle;alignment-baseline:middle">${point.nucleotide}</text>`).join('')
  const labelRects = Array.from(svg.querySelectorAll('text.numbering-label')).map((text) => { const value = text.textContent?.trim() || ''; const x = Number.parseFloat(text.getAttribute('x') || ''); const y = Number.parseFloat(text.getAttribute('y') || ''); const width = Math.max(7, value.length * 4.4 + 2); return Number.isFinite(x) && Number.isFinite(y) ? { left: x - width / 2, right: x + width / 2, top: y - 3.5, bottom: y + 3.5 } : null }).filter((value): value is { left: number; right: number; top: number; bottom: number } => Boolean(value))
  const overlapArea = (left: number, right: number, top: number, bottom: number, rect: { left: number; right: number; top: number; bottom: number }) => Math.max(0, Math.min(right, rect.right) - Math.max(left, rect.left)) * Math.max(0, Math.min(bottom, rect.bottom) - Math.max(top, rect.top))
  const pointToSegmentDistance = (point: { x: number; y: number }, start: { x: number; y: number }, end: { x: number; y: number }) => { const dx = end.x - start.x; const dy = end.y - start.y; const denominator = dx * dx + dy * dy || 1; const ratio = Math.max(0, Math.min(1, ((point.x - start.x) * dx + (point.y - start.y) * dy) / denominator)); return Math.hypot(point.x - (start.x + ratio * dx), point.y - (start.y + ratio * dy)) }
  const avoidanceLines = Array.from(svg.querySelectorAll('line')).map((line) => ({ start: { x: Number.parseFloat(line.getAttribute('x1') || ''), y: Number.parseFloat(line.getAttribute('y1') || '') }, end: { x: Number.parseFloat(line.getAttribute('x2') || ''), y: Number.parseFloat(line.getAttribute('y2') || '') } })).filter((line) => Number.isFinite(line.start.x) && Number.isFinite(line.start.y) && Number.isFinite(line.end.x) && Number.isFinite(line.end.y))
  const orientation = (first: { x: number; y: number }, second: { x: number; y: number }, third: { x: number; y: number }) => (second.x - first.x) * (third.y - first.y) - (second.y - first.y) * (third.x - first.x)
  const segmentsIntersect = (a: { x: number; y: number }, b: { x: number; y: number }, c: { x: number; y: number }, d: { x: number; y: number }) => { const ab = orientation(a, b, c); const ac = orientation(a, b, d); const cd = orientation(c, d, a); const cb = orientation(c, d, b); return ((ab >= 0 && ac <= 0) || (ab <= 0 && ac >= 0)) && ((cd >= 0 && cb <= 0) || (cd <= 0 && cb >= 0)) }
  const lineTouchesRect = (line: { start: { x: number; y: number }; end: { x: number; y: number } }, rect: { left: number; right: number; top: number; bottom: number }) => { const inside = (point: { x: number; y: number }) => point.x >= rect.left && point.x <= rect.right && point.y >= rect.top && point.y <= rect.bottom; if (inside(line.start) || inside(line.end)) return true; const corners = [{ x: rect.left, y: rect.top }, { x: rect.right, y: rect.top }, { x: rect.right, y: rect.bottom }, { x: rect.left, y: rect.bottom }]; return corners.some((corner, index) => segmentsIntersect(line.start, line.end, corner, corners[(index + 1) % corners.length])) }
  const rectDistanceToPoint = (rect: { left: number; right: number; top: number; bottom: number }, point: { x: number; y: number }) => Math.hypot(point.x - Math.max(rect.left, Math.min(point.x, rect.right)), point.y - Math.max(rect.top, Math.min(point.y, rect.bottom)))
  const rectDistanceToLine = (rect: { left: number; right: number; top: number; bottom: number }, line: { start: { x: number; y: number }; end: { x: number; y: number } }) => { if (lineTouchesRect(line, rect)) return 0; const samples = [{ x: rect.left, y: rect.top }, { x: rect.right, y: rect.top }, { x: rect.right, y: rect.bottom }, { x: rect.left, y: rect.bottom }, { x: (rect.left + rect.right) / 2, y: (rect.top + rect.bottom) / 2 }]; return Math.min(...samples.map((sample) => pointToSegmentDistance(sample, line.start, line.end))) }
  const positionLabels = renderCoordinates.filter((point) => point.position % 10 === 0 && !nativeNumberLabels.has(point.position)).map((point) => {
    const vector = { x: point.x - center.x, y: point.y - center.y }; const angle = Math.atan2(vector.y, vector.x); const directions = Array.from({ length: 16 }, (_, index) => { const offset = index === 0 ? 0 : (index % 2 ? 1 : -1) * Math.ceil(index / 2) * Math.PI / 8; return { x: Math.cos(angle + offset), y: Math.sin(angle + offset), priority: Math.abs(offset) } }); const width = Math.max(7, String(point.position).length * 4.4 + 2); let best: any = null
    directions.forEach((direction) => [22, 28, 35, 43, 52].forEach((distance) => { const label = { x: point.x + direction.x * distance, y: point.y + direction.y * distance }; const rect = { left: label.x - width / 2, right: label.x + width / 2, top: label.y - 3.5, bottom: label.y + 3.5 }; const lineStart = { x: point.x + direction.x * 7.6, y: point.y + direction.y * 7.6 }; const lineEnd = { x: point.x + direction.x * (distance - 5.5), y: point.y + direction.y * (distance - 5.5) }; const labelCollisions = labelRects.reduce((total, existing) => total + overlapArea(rect.left, rect.right, rect.top, rect.bottom, existing), 0); const baseClearance = Math.min(...renderCoordinates.filter((candidate) => candidate.position !== point.position).map((candidate) => rectDistanceToPoint(rect, candidate) - 6.6)); const lineCollisions = avoidanceLines.reduce((total, line) => pointToSegmentDistance(line.start, lineStart, lineEnd) < 3.2 || pointToSegmentDistance(line.end, lineStart, lineEnd) < 3.2 ? total + 1 : total, 0); const lineClearance = Math.min(...avoidanceLines.map((line) => rectDistanceToLine(rect, line))); const clearance = Math.min(baseClearance, lineClearance); const unsafe = clearance < 3.5 || labelCollisions > 0 || lineCollisions > 0; const score = unsafe ? -1000000000 + clearance * 100000 - labelCollisions * 10000 - lineCollisions * 1000 - distance : -distance * 100 - direction.priority; if (!best || score > best.score) best = { score, label, rect, lineStart, lineEnd } }))
    labelRects.push(best.rect); avoidanceLines.push({ start: best.lineStart, end: best.lineEnd }); return `<g><line x1="${best.lineStart.x.toFixed(3)}" y1="${best.lineStart.y.toFixed(3)}" x2="${best.lineEnd.x.toFixed(3)}" y2="${best.lineEnd.y.toFixed(3)}" class="numbering-line sequential"/><text class="numbering-label sequential" x="${best.label.x.toFixed(3)}" y="${best.label.y.toFixed(3)}" fill="#000000" font-weight="400" text-anchor="middle" dominant-baseline="middle" alignment-baseline="middle" style="font-weight:400;fill:#000000;text-anchor:middle;dominant-baseline:middle;alignment-baseline:middle">${point.position}</text></g>`
  }).join('')
  const hasTerminalLabels = Array.from(svg.querySelectorAll('text')).some((text) => /^(5|3)'$/.test(text.textContent?.trim() || ''))
  const terminalLabels = !hasTerminalLabels && !hasTitledCoordinates && renderCoordinates.length > 1
    ? (() => {
      const first = renderCoordinates[0]
      const last = renderCoordinates[renderCoordinates.length - 1]
      const label = (point: { x: number; y: number }, value: string) => {
        // R2DT's template-free output uses fixed terminal directions in this
        // layout: 5' below the first base and 3' to the right of the last.
        const direction = value === "5'" ? { x: 0, y: 1 } : { x: 1, y: 0 }
        const distance = 25
        const lineStart = { x: point.x + direction.x * 5.833, y: point.y + direction.y * 5.833 }
        const lineEnd = { x: point.x + direction.x * 19.167, y: point.y + direction.y * 19.167 }
        return `<line class="terminal-line" x1="${lineStart.x.toFixed(3)}" y1="${lineStart.y.toFixed(3)}" x2="${lineEnd.x.toFixed(3)}" y2="${lineEnd.y.toFixed(3)}"/><text class="terminal-label" x="${(point.x + direction.x * distance).toFixed(3)}" y="${(point.y + direction.y * distance).toFixed(3)}" fill="#000000" font-family="Helvetica, Arial, sans-serif" font-size="8.75" font-weight="400" text-anchor="middle" dominant-baseline="middle">${value}</text>`
      }
      return `${label(first, "5'")}${label(last, "3'")}`
    })()
    : ''
  coordinates.forEach((point) => { const group = Array.from(svg.querySelectorAll('g')).find((candidate) => (candidate.querySelector(':scope > title')?.textContent || '').startsWith(`${point.position} (position.label`)); group?.querySelector(':scope > text')?.setAttribute('display', 'none') })
  nucleotideTexts.forEach((text) => text.setAttribute('display', 'none'))
  svg.querySelectorAll('text').forEach((text) => { const value = text.textContent?.trim() || ''; if (text.classList.contains('numbering-label')) { text.setAttribute('display', 'inline'); text.setAttribute('fill', '#000000'); text.setAttribute('font-weight', '400'); text.style.fill = '#000000'; text.style.fontWeight = '400' } if (/^(5|3)'$/.test(value)) { text.setAttribute('fill', '#000000'); text.setAttribute('font-weight', '400'); text.style.fill = '#000000'; text.style.fontWeight = '400' } })
  const viewCoordinates = renderCoordinates
  const structureBounds = viewCoordinates.reduce((bounds, point) => ({ minX: Math.min(bounds.minX, point.x), maxX: Math.max(bounds.maxX, point.x), minY: Math.min(bounds.minY, point.y), maxY: Math.max(bounds.maxY, point.y) }), { minX: Infinity, maxX: -Infinity, minY: Infinity, maxY: -Infinity }); const padding = 58
  if (hasTitledCoordinates || nucleotideTexts.length) {
    svg.setAttribute('viewBox', `${(structureBounds.minX - padding).toFixed(3)} ${(structureBounds.minY - padding).toFixed(3)} ${(structureBounds.maxX - structureBounds.minX + padding * 2).toFixed(3)} ${(structureBounds.maxY - structureBounds.minY + padding * 2).toFixed(3)}`)
  }
  svg.setAttribute('preserveAspectRatio', 'xMidYMid meet')
  const legendWidth = Math.min(240, Math.max(150, structureBounds.maxX - structureBounds.minX - 16))
  const legendRows = Math.ceil(legendItems.length / Math.max(1, Math.floor(legendWidth / 62)))
  const legendHeight = legendRows * 15 + 8
  const structureWidth = structureBounds.maxX - structureBounds.minX + padding * 2
  const legendMarkup = buildAnnotationLegendMarkup(legendItems, structureBounds.minX - padding + (structureWidth - legendWidth) / 2, structureBounds.maxY + padding - legendHeight + 12, legendWidth, 'center')
  renderParent.innerHTML = `${renderParent.innerHTML}<style>.basepairs{stroke:#000000!important;stroke-width:1.09375!important}.backbone{stroke:#cccccc!important;stroke-width:1.09375!important}.numbering-line{fill:none!important;stroke:#cccccc!important;stroke-width:0.546875!important}.terminal-line{fill:none!important;stroke:#cccccc!important;stroke-width:1.09375!important}.numbering-label,.rnameta-marked-base-label,.terminal-label{fill:#000000!important;font-family:Helvetica,Arial,sans-serif!important;font-size:8.75px!important;font-weight:400!important}</style><g class="rnameta-terminal-labels" pointer-events="none">${terminalLabels}</g><g class="rnameta-base-nodes" pointer-events="none">${baseNodes}</g><g class="rnameta-marked-base-labels" pointer-events="none">${markedLetters}</g><g class="rnameta-position-labels" pointer-events="none">${positionLabels}</g>${legendMarkup}`
  return svg.outerHTML
}

export function StructureViewer({ src, markup, highlights = [] }: { src: string; markup?: string; highlights?: StructureHighlight[] }) {
  const viewportRef = useRef<HTMLDivElement>(null); const dragRef = useRef<{ x: number; y: number; offset: { x: number; y: number } } | null>(null); const [scale, setScale] = useState(1); const [offset, setOffset] = useState({ x: 0, y: 0 })
  useEffect(() => { setScale(1); setOffset({ x: 0, y: 0 }) }, [markup, src])
  useEffect(() => { const viewport = viewportRef.current; if (!viewport) return undefined; const handleWheel = (event: WheelEvent) => { event.preventDefault(); event.stopPropagation(); const rect = viewport.getBoundingClientRect(); const anchor = { x: event.clientX - rect.left, y: event.clientY - rect.top }; setScale((current) => { const next = Math.min(6, Math.max(.05, current * (event.deltaY < 0 ? 1.12 : .89))); if (next === current) return current; const ratio = next / current; const delta = { x: anchor.x - rect.width / 2, y: anchor.y - rect.height / 2 }; setOffset((previous) => ({ x: ratio * previous.x + (1 - ratio) * delta.x, y: ratio * previous.y + (1 - ratio) * delta.y })); return next }) }; viewport.addEventListener('wheel', handleWheel, { passive: false }); return () => viewport.removeEventListener('wheel', handleWheel) }, [])
  const renderedMarkup = useMemo(() => buildStructureMarkup(markup || '', highlights), [markup, highlights])
  const handlePointerDown = (event: React.PointerEvent<HTMLDivElement>) => { event.preventDefault(); event.currentTarget.setPointerCapture(event.pointerId); dragRef.current = { x: event.clientX, y: event.clientY, offset } }
  const handlePointerMove = (event: React.PointerEvent<HTMLDivElement>) => { if (!dragRef.current) return; setOffset({ x: dragRef.current.offset.x + event.clientX - dragRef.current.x, y: dragRef.current.offset.y + event.clientY - dragRef.current.y }) }
  return <div className="rna-annotation-structure-viewer"><div ref={viewportRef} className="rna-annotation-structure-viewport" onPointerDown={handlePointerDown} onPointerMove={handlePointerMove} onPointerUp={() => { dragRef.current = null }} onPointerCancel={() => { dragRef.current = null }}>{renderedMarkup ? <div className="rna-annotation-structure-inline" style={{ transform: `translate(${offset.x}px, ${offset.y}px) scale(${scale})` }} dangerouslySetInnerHTML={{ __html: renderedMarkup }} /> : <img className="rna-annotation-structure-preview__image" src={src} alt="R2DT RNA secondary structure" style={{ transform: `translate(${offset.x}px, ${offset.y}px) scale(${scale})` }} draggable="false" />}</div></div>
}

export function AnnotationTracks({ highlights = [], sequenceLength = 0 }: { highlights?: StructureHighlight[]; sequenceLength?: number }) {
  const tracks = useMemo(() => { const labels = [...new Set(highlights.map((annotation) => String(annotation.label || annotation.category || 'annotation').trim()))].sort((a, b) => a.localeCompare(b)); const usedColors = new Set<string>(); return labels.map((label) => { const key = label.toLowerCase(); const color = colorForLabel(key, usedColors); usedColors.add(color); const lanes: Array<Array<{ start: number; end: number; display: string; labelLeft: number; labelRight: number }>> = []; highlights.filter((annotation) => String(annotation.label || annotation.category || 'annotation').trim() === label).map((annotation) => { const start = Number(annotation.start); const end = Number(annotation.end); const display = start === end ? `${label}(${start})` : `${label}(${start}-${end})`; const x = ((start - 1) / Math.max(1, sequenceLength - 1)) * 816; const xEnd = ((end - 1) / Math.max(1, sequenceLength - 1)) * 816; const textWidth = Math.max(36, display.length * 5.1); return { start, end, display, labelLeft: (x + xEnd) / 2 - textWidth / 2, labelRight: (x + xEnd) / 2 + textWidth / 2 } }).sort((a, b) => a.start - b.start || a.end - b.end).forEach((interval) => { const laneIndex = lanes.findIndex((lane) => lane[lane.length - 1].labelRight + 5 < interval.labelLeft); if (laneIndex === -1) lanes.push([interval]); else lanes[laneIndex].push(interval) }); return { label, color, lanes } }) }, [highlights, sequenceLength])
  if (!tracks.length || !sequenceLength) return null
  const left = 24; const right = 24; const width = 1000; const laneHeight = 27; const trackGap = 18; const top = 45; const rows = tracks.reduce((count, track) => count + track.lanes.length, 0); const height = top + rows * laneHeight + Math.max(0, tracks.length - 1) * trackGap + 31; const axisY = height - 21; const tickPositions = Array.from({ length: Math.floor(sequenceLength / 10) }, (_, index) => (index + 1) * 10); let y = top; const x = (position: number) => left + ((position - 1) / Math.max(1, sequenceLength - 1)) * (width - left - right)
  return <div className="rna-annotation-tracks"><svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label="RNA annotation tracks"><g className="rna-annotation-tracks__legend" aria-label="Annotation color legend">{tracks.map((track, index) => <g key={`legend-${track.label}`} transform={`translate(${left + index * 110}, 16)`}><circle cx="0" cy="0" r="3.3" fill={track.color} stroke={darkerColor(track.color)} strokeWidth="0.8" /><text x="6" y="3.5" className="rna-annotation-tracks__legend-label">{track.label}</text></g>)}</g>{tracks.map((track) => { const trackTop = y; const trackHeight = track.lanes.length * laneHeight; y += trackHeight + trackGap; return <g key={track.label}><rect className="rna-annotation-tracks__band" x={left - 18} y={trackTop - 7} width={width - left - right + 36} height={trackHeight + 14} rx="4" />{tickPositions.map((position) => <line className="rna-annotation-tracks__grid" key={`${track.label}-${position}`} x1={x(position)} x2={x(position)} y1={trackTop - 3} y2={trackTop + trackHeight + 3} />)}{track.lanes.map((lane, laneIndex) => { const laneY = trackTop + laneIndex * laneHeight; return <g key={`${track.label}-${laneIndex}`}><line className="rna-annotation-tracks__baseline" x1={left} x2={width - right} y1={laneY + 19} y2={laneY + 19} />{lane.map((interval, intervalIndex) => { const start = x(interval.start); const end = x(interval.end); return <g key={`${interval.start}-${interval.end}-${intervalIndex}`}><text className="rna-annotation-tracks__interval" x={(start + end) / 2} y={laneY + 7} textAnchor="middle">{interval.display}</text><rect x={start} y={laneY + 16} width={Math.max(5, end - start + 2)} height="6" rx="2" fill={track.color} stroke={darkerColor(track.color)} strokeWidth="1" /></g> })}</g> })}</g> })}<line className="rna-annotation-tracks__axis" x1={left} x2={width - right} y1={axisY} y2={axisY} />{tickPositions.map((position) => <g key={position}><line className="rna-annotation-tracks__tick" x1={x(position)} x2={x(position)} y1={axisY} y2={axisY + 5} /><text className="rna-annotation-tracks__axis-label" x={x(position)} y={height - 6} textAnchor="middle">{position}</text></g>)}<text className="rna-annotation-tracks__axis-label" x={left} y={height - 6} textAnchor="start">1</text><text className="rna-annotation-tracks__axis-label" x={width - right} y={height - 6} textAnchor="end">{sequenceLength}</text></svg></div>
}
