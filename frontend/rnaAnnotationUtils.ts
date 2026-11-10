export type StructureHighlight = { start?: number; end?: number; label?: string; category?: string }

function darkerColor(color: string) {
  const match = String(color || '').match(/^#([0-9a-f]{6})$/i)
  if (!match) return '#50504c'
  const values = [0, 2, 4].map((offset) => Number.parseInt(match[1].slice(offset, offset + 2), 16))
  return `#${values.map((value) => Math.round(value * .64).toString(16).padStart(2, '0')).join('')}`
}

function escapeXml(value: string) {
  return value.replace(/[&<>"']/g, (character) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&apos;' }[character] || character))
}

export function buildAnnotationLegendMarkup(items: Array<{ label: string; color: string }>, x: number, y: number, maxWidth = 180, align: 'start' | 'center' = 'start') {
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
