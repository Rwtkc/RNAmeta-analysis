meta_plot_controls_config <- function(id) {
  ns <- NS(id)

  list(
    ids = list(
      pltTxType = ns("pltTxType"),
      headOrtail = ns("headOrtail"),
      txPrimaryOnly = ns("txPrimaryOnly"),
      txlncrnaOverlapmrna = ns("txlncrnaOverlapmrna"),
      mapFilterTranscript = ns("mapFilterTranscript"),
      enableCI = ns("enableCI"),
      stSampleModle = ns("stSampleModle"),
      stSampleNum = ns("stSampleNum"),
      stAmblguity = ns("stAmblguity"),
      txfiveutrMinLength = ns("txfiveutrMinLength"),
      txcdsMinLength = ns("txcdsMinLength"),
      txthreeutrMinLength = ns("txthreeutrMinLength"),
      txlongNcrnaMinLength = ns("txlongNcrnaMinLength"),
      txpromoterLength = ns("txpromoterLength"),
      txtailLength = ns("txtailLength"),
      adjust = ns("adjust"),
      overlapIndex = ns("overlapIndex"),
      siteLengthIndex = ns("siteLengthIndex")
    ),
    defaults = list(
      pltTxType = "tx",
      headOrtail = "TRUE",
      txPrimaryOnly = "TRUE",
      txlncrnaOverlapmrna = "TRUE",
      mapFilterTranscript = "TRUE",
      enableCI = "FALSE",
      stSampleModle = "Equidistance",
      stSampleNum = 10,
      stAmblguity = 5,
      txfiveutrMinLength = 100,
      txcdsMinLength = 100,
      txthreeutrMinLength = 100,
      txlongNcrnaMinLength = 100,
      txpromoterLength = 1000,
      txtailLength = 1000,
      adjust = 1,
      overlapIndex = 1,
      siteLengthIndex = 1
    ),
    choices = list(
      transcriptType = list(
        list(label = "exon", value = "tx"),
        list(label = "mRNA", value = "mrna"),
        list(label = "lncRNA", value = "ncrna")
      ),
      samplingMethods = list(
        list(label = "\"Equidistance\"", value = "Equidistance"),
        list(label = "\"random\"", value = "random")
      ),
      binary = list(
        list(label = "TRUE", value = "TRUE"),
        list(label = "FALSE", value = "FALSE")
      )
    )
  )
}

meta_plot_control_values <- function(input, defaults) {
  list(
    pltTxType = if (is.null(input$pltTxType) || identical(input$pltTxType, "")) defaults$pltTxType else input$pltTxType,
    headOrtail = if (is.null(input$headOrtail) || identical(input$headOrtail, "")) defaults$headOrtail else input$headOrtail,
    txPrimaryOnly = if (is.null(input$txPrimaryOnly) || identical(input$txPrimaryOnly, "")) defaults$txPrimaryOnly else input$txPrimaryOnly,
    txlncrnaOverlapmrna = if (is.null(input$txlncrnaOverlapmrna) || identical(input$txlncrnaOverlapmrna, "")) defaults$txlncrnaOverlapmrna else input$txlncrnaOverlapmrna,
    mapFilterTranscript = if (is.null(input$mapFilterTranscript) || identical(input$mapFilterTranscript, "")) defaults$mapFilterTranscript else input$mapFilterTranscript,
    enableCI = if (is.null(input$enableCI) || identical(input$enableCI, "")) defaults$enableCI else input$enableCI,
    stSampleModle = if (is.null(input$stSampleModle) || identical(input$stSampleModle, "")) defaults$stSampleModle else input$stSampleModle,
    stSampleNum = if (is.null(input$stSampleNum) || identical(input$stSampleNum, "")) defaults$stSampleNum else suppressWarnings(as.integer(input$stSampleNum)),
    stAmblguity = if (is.null(input$stAmblguity) || identical(input$stAmblguity, "")) defaults$stAmblguity else suppressWarnings(as.integer(input$stAmblguity)),
    txfiveutrMinLength = if (is.null(input$txfiveutrMinLength) || identical(input$txfiveutrMinLength, "")) defaults$txfiveutrMinLength else suppressWarnings(as.integer(input$txfiveutrMinLength)),
    txcdsMinLength = if (is.null(input$txcdsMinLength) || identical(input$txcdsMinLength, "")) defaults$txcdsMinLength else suppressWarnings(as.integer(input$txcdsMinLength)),
    txthreeutrMinLength = if (is.null(input$txthreeutrMinLength) || identical(input$txthreeutrMinLength, "")) defaults$txthreeutrMinLength else suppressWarnings(as.integer(input$txthreeutrMinLength)),
    txlongNcrnaMinLength = if (is.null(input$txlongNcrnaMinLength) || identical(input$txlongNcrnaMinLength, "")) defaults$txlongNcrnaMinLength else suppressWarnings(as.integer(input$txlongNcrnaMinLength)),
    txpromoterLength = if (is.null(input$txpromoterLength) || identical(input$txpromoterLength, "")) defaults$txpromoterLength else suppressWarnings(as.integer(input$txpromoterLength)),
    txtailLength = if (is.null(input$txtailLength) || identical(input$txtailLength, "")) defaults$txtailLength else suppressWarnings(as.integer(input$txtailLength)),
    adjust = if (is.null(input$adjust) || identical(input$adjust, "")) defaults$adjust else suppressWarnings(as.numeric(input$adjust)),
    overlapIndex = if (is.null(input$overlapIndex) || identical(input$overlapIndex, "")) defaults$overlapIndex else suppressWarnings(as.numeric(input$overlapIndex)),
    siteLengthIndex = if (is.null(input$siteLengthIndex) || identical(input$siteLengthIndex, "")) defaults$siteLengthIndex else suppressWarnings(as.numeric(input$siteLengthIndex))
  )
}
