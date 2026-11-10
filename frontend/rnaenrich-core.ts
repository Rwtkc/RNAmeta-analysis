export type MatchMode = 'exact' | 'interval'
export type BackgroundMode = 'none' | 'union' | 'upload'

export type BedRecord = {
  chrom: string
  start: number
  end: number
  sourceLine: string
}

export type ParsedBed = {
  name: string
  records: BedRecord[]
  rawCount: number
  duplicateCount: number
  invalidLines: { line: number; content: string; reason: string }[]
}

export type PairResult = {
  left: string
  right: string
  leftSize: number
  rightSize: number
  overlapCount: number
  leftMatched: number
  rightMatched: number
  leftRate: number
  rightRate: number
  jaccard: number
  overlapBp: number
  unionBp: number
  pValue: number | null
  fdr: number | null
}

const canonicalChrom = (value: string) => {
  const trimmed = value.trim()
  if (!trimmed) return trimmed
  const withoutPrefix = trimmed.replace(/^chr/i, '')
  if (/^(\d+|x|y|m|mt)$/i.test(withoutPrefix)) {
    const normalized = withoutPrefix.toUpperCase() === 'MT' ? 'M' : withoutPrefix.toUpperCase()
    return `chr${normalized}`
  }
  return trimmed
}

export const recordKey = (record: BedRecord) => `${record.chrom}\t${record.start}\t${record.end}`

export function parseBed(name: string, text: string): ParsedBed {
  const records: BedRecord[] = []
  const invalidLines: ParsedBed['invalidLines'] = []
  const seen = new Set<string>()
  let rawCount = 0
  let duplicateCount = 0

  text.split(/\r?\n/).forEach((sourceLine, index) => {
    const line = sourceLine.trim()
    if (!line || line.startsWith('#') || line.startsWith('track ') || line.startsWith('browser ')) return
    rawCount += 1
    const fields = line.split(/\s+/)
    const start = Number(fields[1])
    const end = Number(fields[2])
    if (fields.length < 3) {
      invalidLines.push({ line: index + 1, content: sourceLine, reason: 'Expected at least three columns.' })
      return
    }
    if (!Number.isSafeInteger(start) || !Number.isSafeInteger(end) || start < 0 || end <= start) {
      invalidLines.push({ line: index + 1, content: sourceLine, reason: 'Coordinates must be integers with 0 ≤ start < end.' })
      return
    }
    const record = { chrom: canonicalChrom(fields[0]), start, end, sourceLine }
    const key = recordKey(record)
    if (seen.has(key)) {
      duplicateCount += 1
      return
    }
    seen.add(key)
    records.push(record)
  })

  records.sort((a, b) => a.chrom.localeCompare(b.chrom, undefined, { numeric: true }) || a.start - b.start || a.end - b.end)
  return { name, records, rawCount, duplicateCount, invalidLines }
}

function groupByChrom(records: BedRecord[]) {
  const grouped = new Map<string, BedRecord[]>()
  records.forEach((record) => {
    const current = grouped.get(record.chrom) ?? []
    current.push(record)
    grouped.set(record.chrom, current)
  })
  return grouped
}

function intervalHits(left: BedRecord[], right: BedRecord[]) {
  const leftByChrom = groupByChrom(left)
  const rightByChrom = groupByChrom(right)
  let overlapCount = 0
  let leftMatched = 0
  let rightMatched = 0

  for (const [chrom, leftRecords] of leftByChrom) {
    const rightRecords = rightByChrom.get(chrom) ?? []
    const matchedRight = new Set<number>()
    let firstCandidate = 0
    leftRecords.forEach((a) => {
      while (firstCandidate < rightRecords.length && rightRecords[firstCandidate].end <= a.start) firstCandidate += 1
      let matched = false
      for (let rightIndex = firstCandidate; rightIndex < rightRecords.length; rightIndex += 1) {
        const b = rightRecords[rightIndex]
        if (b.start >= a.end) break
        if (b.end > a.start) {
          overlapCount += 1
          matched = true
          matchedRight.add(rightIndex)
        }
      }
      if (matched) leftMatched += 1
    })
    rightMatched += matchedRight.size
  }
  return { overlapCount, leftMatched, rightMatched }
}

function mergeRecords(records: BedRecord[]) {
  const grouped = groupByChrom(records)
  const merged: BedRecord[] = []
  for (const [chrom, chromRecords] of grouped) {
    let currentStart = -1
    let currentEnd = -1
    chromRecords.forEach((record) => {
      if (currentStart < 0 || record.start > currentEnd) {
        if (currentStart >= 0) merged.push({ chrom, start: currentStart, end: currentEnd, sourceLine: '' })
        currentStart = record.start
        currentEnd = record.end
      } else {
        currentEnd = Math.max(currentEnd, record.end)
      }
    })
    if (currentStart >= 0) merged.push({ chrom, start: currentStart, end: currentEnd, sourceLine: '' })
  }
  return merged
}

function mergedLength(records: BedRecord[]) {
  return mergeRecords(records).reduce((total, record) => total + record.end - record.start, 0)
}

function intersectionLength(left: BedRecord[], right: BedRecord[]) {
  const leftByChrom = groupByChrom(mergeRecords(left))
  const rightByChrom = groupByChrom(mergeRecords(right))
  let total = 0
  for (const [chrom, leftRecords] of leftByChrom) {
    const rightRecords = rightByChrom.get(chrom)
    if (!rightRecords) continue
    let i = 0
    let j = 0
    while (i < leftRecords.length && j < rightRecords.length) {
      const a = leftRecords[i]
      const b = rightRecords[j]
      const start = Math.max(a.start, b.start)
      const end = Math.min(a.end, b.end)
      if (start < end) total += end - start
      if (a.end <= b.end) i += 1
      else j += 1
    }
  }
  return total
}

function logGamma(value: number): number {
  const coefficients = [676.5203681218851, -1259.1392167224028, 771.3234287776531, -176.6150291621406, 12.507343278686905, -0.13857109526572012, 9.984369578019572e-6, 1.5056327351493116e-7]
  if (value < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * value)) - logGamma(1 - value)
  const z = value - 1
  let x = 0.9999999999998099
  coefficients.forEach((coefficient, index) => { x += coefficient / (z + index + 1) })
  const t = z + coefficients.length - 0.5
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x)
}

function logChoose(n: number, k: number) {
  if (k < 0 || k > n) return Number.NEGATIVE_INFINITY
  return logGamma(n + 1) - logGamma(k + 1) - logGamma(n - k + 1)
}

export function hypergeometricUpperTail(population: number, successStates: number, draws: number, observed: number) {
  const maximum = Math.min(successStates, draws)
  const minimum = Math.max(0, draws - (population - successStates))
  if (observed <= minimum) return 1
  if (observed > maximum) return 0
  const denominator = logChoose(population, draws)
  let probability = Math.exp(logChoose(successStates, observed) + logChoose(population - successStates, draws - observed) - denominator)
  let total = probability
  for (let value = observed; value < maximum; value += 1) {
    const numerator = (successStates - value) * (draws - value)
    const divisor = (value + 1) * (population - successStates - draws + value + 1)
    if (divisor === 0) break
    probability *= numerator / divisor
    total += probability
  }
  return Math.min(1, total)
}

export function adjustBenjaminiHochberg(values: number[]) {
  const indexed = values.map((value, index) => ({ value, index })).sort((a, b) => a.value - b.value)
  const adjusted = Array(values.length).fill(1) as number[]
  let previous = 1
  for (let rank = indexed.length; rank >= 1; rank -= 1) {
    const item = indexed[rank - 1]
    previous = Math.min(previous, item.value * indexed.length / rank)
    adjusted[item.index] = Math.min(1, previous)
  }
  return adjusted
}

export function compareBeds(files: ParsedBed[], mode: MatchMode, backgroundSize: number | null): PairResult[] {
  const pairs: PairResult[] = []
  for (let i = 0; i < files.length; i += 1) {
    for (let j = i + 1; j < files.length; j += 1) {
      const left = files[i]
      const right = files[j]
      if (mode === 'exact') {
        const rightKeys = new Set(right.records.map(recordKey))
        const overlapCount = left.records.reduce((count, record) => count + Number(rightKeys.has(recordKey(record))), 0)
        const unionCount = left.records.length + right.records.length - overlapCount
        const pValue = backgroundSize === null ? null : hypergeometricUpperTail(backgroundSize, right.records.length, left.records.length, overlapCount)
        pairs.push({
          left: left.name, right: right.name, leftSize: left.records.length, rightSize: right.records.length,
          overlapCount, leftMatched: overlapCount, rightMatched: overlapCount,
          leftRate: left.records.length ? overlapCount / left.records.length : 0,
          rightRate: right.records.length ? overlapCount / right.records.length : 0,
          jaccard: unionCount ? overlapCount / unionCount : 0,
          overlapBp: overlapCount, unionBp: unionCount, pValue, fdr: null,
        })
      } else {
        const hits = intervalHits(left.records, right.records)
        const leftBp = mergedLength(left.records)
        const rightBp = mergedLength(right.records)
        const overlapBp = intersectionLength(left.records, right.records)
        const unionBp = leftBp + rightBp - overlapBp
        pairs.push({
          left: left.name, right: right.name, leftSize: left.records.length, rightSize: right.records.length,
          ...hits,
          leftRate: left.records.length ? hits.leftMatched / left.records.length : 0,
          rightRate: right.records.length ? hits.rightMatched / right.records.length : 0,
          jaccard: unionBp ? overlapBp / unionBp : 0,
          overlapBp, unionBp, pValue: null, fdr: null,
        })
      }
    }
  }
  const pValues = pairs.flatMap((pair) => pair.pValue === null ? [] : [pair.pValue])
  if (pValues.length) {
    const adjusted = adjustBenjaminiHochberg(pValues)
    let index = 0
    pairs.forEach((pair) => {
      if (pair.pValue !== null) {
        pair.fdr = adjusted[index]
        index += 1
      }
    })
  }
  return pairs
}

export function exactMembership(files: ParsedBed[]) {
  const memberships = new Map<string, Set<string>>()
  files.forEach((file) => file.records.forEach((record) => {
    const key = recordKey(record)
    const names = memberships.get(key) ?? new Set<string>()
    names.add(file.name)
    memberships.set(key, names)
  }))
  const combinations = new Map<string, { count: number; members: string[] }>()
  memberships.forEach((names) => {
    const members = files.map((file) => file.name).filter((name) => names.has(name))
    const key = JSON.stringify(members)
    const current = combinations.get(key)
    combinations.set(key, { count: (current?.count ?? 0) + 1, members })
  })
  return [...combinations.values()]
    .map(({ count, members }) => ({ label: members.join(' ∩ '), count, members }))
    .sort((a, b) => b.count - a.count || a.label.localeCompare(b.label))
}

export function exactSharedRecords(left: ParsedBed, right: ParsedBed) {
  const rightKeys = new Set(right.records.map(recordKey))
  return left.records.filter((record) => rightKeys.has(recordKey(record)))
}
