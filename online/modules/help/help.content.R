help_sections <- function() {
  list(
    list(
      id = "overview",
      label = "Overview",
      summary = paste(
        "Understand what RNAmeta Online covers, how it relates to the desktop",
        "client, and which analysis surfaces are currently available online."
      ),
      blocks = list(
        list(
          title = "What RNAmeta Online is designed for",
          paragraphs = c(
            paste(
              "RNAmeta Online is a browser-based workspace for interpreting BED-based",
              "RNA modification peak data. It retains the main analysis modules",
              "used across RNAmeta Online so that transcript-scale, gene-level, exon-level,",
              "site-level, and motif-level results can be examined within one workflow."
            ),
            paste(
              "RNAmeta Online is organized around a saved Upload / Run context.",
              "Once you save species information and one or more BED files, the",
              "downstream modules reuse the same dataset so that results can be",
              "compared across views without repeated uploading."
            )
          ),
          lists = list(
            list(
              label = "Current web analysis coverage",
              items = c(
                "Upload / Run for species context, BED saving, preview, and demo loading",
                "Peak Distribution for promoter, UTR, CDS, intron, and intergenic annotation summaries",
                "Genome Browser for IGV-Web inspection of saved BED intervals against reference tracks",
                  "Structure for transcript mapping and RNAfold-based secondary structure jobs",
                "Gene Statistics for Gene Type, Peak Gene Size, and Gene Matrix",
                "Exon Statistics for Peak Exon Size, Peak Exon Type, and Peak Exon Num",
                "Site for Transcription, Translation, and Splicesite boundary analyses",
                "Motif for background HOMER, MEME, and STREME jobs"
              )
            )
          )
        ),
        list(
          title = "How the web app relates to RNAmeta Desktop",
          paragraphs = c(
            paste(
              "The desktop client and the public web app share the same overall",
              "analysis vocabulary, so Overview, Upload / Run, Peak Distribution,",
              "Genome Browser, Structure, Gene Statistics, Exon Statistics, Site,",
              "Export, and Troubleshooting remain natural reference points on both sides."
            ),
            paste(
              "They are not identical products. The desktop client includes a",
              "local-client-oriented Help system and desktop runtime behavior,",
              "while the web app starts directly from Upload / Run and also",
              "exposes a Motif job workflow that is not part of the current",
              "desktop navigation."
            )
          ),
          lists = list(
            list(
              label = "Practical scope notes",
              items = c(
                "Use this Help page for the current web workflow rather than for local desktop setup.",
                "Treat Motif as a web-specific module in the current RNAmeta Online release.",
                "When the saved Upload / Run context changes, rerun the affected downstream modules."
              )
            )
          )
        )
      )
    ),
    list(
      id = "getting-started",
      label = "Getting Started",
      summary = paste(
        "Follow the saved-session workflow from input preparation to module-specific",
        "analysis and result interpretation."
      ),
      blocks = list(
        list(
          title = "Step 1: Save an input session",
          paragraphs = c(
            paste(
              "Begin in Upload / Run and save a BED session. This session provides",
              "the shared input context for all downstream analysis pages."
            )
          )
        ),
        list(
          title = "Step 2: Choose the analysis that answers your question",
          paragraphs = c(
            paste(
              "Use Peak Distribution for feature annotation, Genome Browser for",
              "locus inspection, Structure,",
              "and the statistics, site, or motif pages for their respective analyses."
            )
          )
        ),
        list(
          title = "Step 3: Rerun after changing the saved session",
          paragraphs = c(
            paste(
              "When you save a different BED context or switch species, treat the",
              "previous results as outdated and rerun the modules before",
              "interpreting or exporting their outputs."
            )
          )
        )
      )
    ),
    list(
      id = "upload-run",
      label = "Upload / Run",
      summary = paste(
        "Save the species-aware BED session used by every downstream page and",
        "inspect the uploaded file preview before analysis."
      ),
      blocks = list(
        list(
          title = "What Upload / Run controls",
          paragraphs = c(
            paste(
              "Upload / Run is the entry point for the shared web session. It",
              "captures the current data source choice, species label, and saved",
              "BED files so the rest of the app can work from one consistent",
              "context."
            ),
            paste(
              "Although the form includes a data-source selector, the current",
              "downstream workflow is centered on saved BED files. A valid saved BED",
              "context is what enables the analysis pages."
            )
          )
        ),
        list(
          title = "Working with one or more BED files",
          paragraphs = c(
            paste(
              "RNAmeta Online accepts BED-formatted chromosomal coordinates from common",
              "peak-calling workflows. Saving multiple BED files is useful when you",
              "want cross-sample comparisons in modules such as Peak Distribution,",
              "Gene Statistics, Exon Statistics, Site, or Motif."
            ),
            paste(
              "File labels remain visible in legends, tooltips, sample selectors,",
              "and exports, so clear file naming improves interpretation."
            )
          ),
          lists = list(
            list(
              label = "Recommended habits",
              items = c(
                "Keep one saved session biologically coherent by using the same species context across files.",
                "Use distinct sample names because those labels reappear in charts and export tables.",
                "Resave the session when you switch to a new sample batch instead of interpreting old results."
              )
            )
          )
        ),
        list(
          title = "Preview, demo, and downstream dependencies",
          paragraphs = c(
            paste(
              "After Save, the Session Summary shows the current species, saved-file",
              "count, and a BED preview. This preview is designed as a quick format",
              "check rather than a full browser."
            ),
            paste(
              "Every downstream module depends on this saved context. If Upload /",
              "Run is incomplete, most analysis buttons stay disabled and the page",
              "shows a hint telling you to save BED files first."
            )
          )
        )
      )
    ),
    list(
      id = "peak-distribution",
      label = "Peak Distribution",
      summary = paste(
        "Annotate saved BED intervals against transcript-related feature classes",
        "and compare feature-level composition across files."
      ),
      blocks = list(
        list(
          title = "What Peak Distribution answers",
          paragraphs = c(
            paste(
              "Peak Distribution shows where the current BED intervals land relative",
              "to promoter, UTR, CDS, intron, codon, and intergenic-style feature",
              "classes. It is the page to use when your main question is functional",
              "annotation rather than transcript-relative shape."
            ),
            paste(
              "This module works well near the start of interpretation because it",
              "quickly reveals whether the dataset is dominated by a few feature",
              "classes or spread broadly across gene architecture."
            )
          )
        ),
        list(
          title = "Display settings and export behavior",
          paragraphs = c(
            paste(
              "The sidebar lets you choose which feature classes stay visible in the",
              "current plot. That visible-state filter matters for export as well:",
              "the exported table follows the active rendered view instead of an",
              "invisible superset of categories."
            ),
            paste(
              "If the result looks unusual, first confirm that the Upload / Run",
              "species context matches the BED files you saved, because annotation",
              "mismatch will distort feature assignment."
            )
          )
        )
      )
    ),
    list(
      id = "genome-browser",
      label = "Genome Browser",
      summary = paste(
        "Inspect the saved Upload / Run BED intervals in an IGV-Web workspace",
        "with species-matched reference and annotation tracks."
      ),
      blocks = list(
        list(
          title = "What Genome Browser is useful for",
          paragraphs = c(
            paste(
              "Genome Browser is the locus-level inspection page for RNAmeta Online.",
              "Use it when a global summary suggests an interesting pattern and",
              "you want to confirm how saved BED intervals sit on the genome",
              "together with reference sequence and annotation tracks."
            ),
            paste(
              "The page reuses the saved Upload / Run species and BED context.",
              "After the browser loads, use the IGV-Web navigation controls to",
              "search loci, zoom, pan, and visually compare uploaded intervals",
              "with the selected reference annotation."
            )
          )
        ),
        list(
          title = "Required input and browser refresh",
          paragraphs = c(
            paste(
              "Save at least one BED file in Upload / Run before opening Genome",
              "Browser. The selected species controls which reference files and",
              "annotation tracks are loaded, so a species mismatch will make the",
              "browser view misleading even if the BED file itself is valid."
            ),
            paste(
              "Use Refresh Browser after changing the saved Upload / Run context",
              "or if the browser frame is still showing an older track set. The",
              "snapshot panel above the output summarizes the active species,",
              "saved BED files, and current locus context."
            )
          ),
          lists = list(
            list(
              label = "Coordinate habits",
              items = c(
                "Use chromosome names that match the selected reference, such as chr1-style or species-specific chromosome labels.",
                "Check BED file names and track labels before interpreting overlapping intervals.",
                "Refresh the browser after saving a different BED batch or switching species."
              )
            )
          )
        )
      )
    ),
    list(
      id = "structure",
      label = "Structure",
      summary = paste(
        "Map saved BED intervals to transcripts, submit selected rows to",
        "RNAfold, and reopen completed Structure jobs by job id."
      ),
      blocks = list(
        list(
          title = "What Structure does",
          paragraphs = c(
            paste(
              "Run Structure first to map saved BED intervals onto transcript",
              "coordinates. The resulting table keeps the first three BED columns",
              "and adds Gene ID, Transcript ID, and Transcript interval."
            ),
            paste(
              "Select a mapped row and choose Submit to RNAfold. The server derives",
              "the reference sequence for that transcript interval and renders an",
              "interactive RNA secondary-structure preview in the output panel."
            )
          )
        ),
        list(
          title = "Mapping table and RNAfold submission",
          paragraphs = c(
            paste(
              "Structure does not require manual coordinate entry. After mapping",
              "finishes, use the table search and pagination controls to find a",
              "transcript or BED row of interest."
            ),
            paste(
              "Confirm the selected row in the RNAfold dialog before submitting.",
              "The job receives the mapped BED interval and uses the full transcript",
              "window for the current Structure preview workflow."
            )
          ),
          lists = list(
            list(
              label = "Before clicking Run Structure",
              items = c(
                "Save a BED file in Upload / Run and choose the same species used to generate the BED coordinates.",
                "Run Structure and wait for the Transcript Mapping Table to appear.",
                "Search or paginate to a mapped row, then confirm its RNAfold submission.",
                "Keep the Structure Job ID if you want to query the result after leaving the page."
              )
            )
          )
        ),
        list(
          title = "Background jobs and result lookup",
          paragraphs = c(
            paste(
              "Structure jobs run on the server workspace. The Structure Job panel",
              "shows the shared server queue, current Job ID, and live queued, running,",
              "failed, or completed status. Completed results can be reopened from the",
              "Structure Job ID field."
            ),
            paste(
              "Structure results follow the same short-term web-result retention",
              "model used by server-side jobs in RNAmeta Online. Keep the job id",
              "while you are still reviewing the result, and export figures you",
              "need to preserve."
            )
          )
        )
      )
    ),
    list(
      id = "gene-statistics",
      label = "Gene Statistics",
      summary = paste(
        "Use Gene Type, Peak Gene Size, and Gene Matrix to summarize biotype,",
        "length-linked structure, and cross-sample overlap."
      ),
      blocks = list(
        list(
          title = "Gene Type",
          paragraphs = c(
            paste(
              "Gene Type summarizes the transcript biotypes associated with the",
              "saved peak set. Use it when you want a quick answer to whether the",
              "signals are concentrated in coding transcripts, non-coding classes,",
              "or a mixed background."
            )
          )
        ),
        list(
          title = "Peak Gene Size",
          paragraphs = c(
            paste(
              "Peak Gene Size focuses on transcript-length distributions for the",
              "genes linked to the saved peaks. It is useful when you suspect that",
              "signal enrichment is biased toward shorter or longer transcript",
              "structures."
            ),
            paste(
              "Treat the output as a structural distribution view and export the",
              "table if you want to perform additional statistics elsewhere."
            )
          )
        ),
        list(
          title = "Gene Matrix",
          paragraphs = c(
            paste(
              "Gene Matrix is the multi-file overlap view. It requires 2 to 5 saved",
              "BED files and compares which genes are shared or unique across the",
              "current session."
            ),
            paste(
              "Use it when the main question is concordance between samples rather",
              "than within-sample structure."
            )
          ),
          lists = list(
            list(
              label = "When Gene Matrix is the right first stop",
              items = c(
                "You need to compare replicate-like BED files side by side.",
                "You want a shared-versus-unique summary before running exon or site-level pages.",
                "You have between 2 and 5 saved BED files in the current Upload / Run session."
              )
            )
          )
        )
      )
    ),
    list(
      id = "exon-statistics",
      label = "Exon Statistics",
      summary = paste(
        "Inspect exon-length, exon-position, and exon-count structure for the",
        "transcripts associated with the saved peaks."
      ),
      blocks = list(
        list(
          title = "Peak Exon Size",
          paragraphs = c(
            paste(
              "Peak Exon Size summarizes exon-length distributions while preserving",
              "the distinction between first, middle, and last exon classes. Use",
              "it when exon length might explain the structural pattern behind your",
              "peak-associated transcripts."
            )
          )
        ),
        list(
          title = "Peak Exon Type",
          paragraphs = c(
            paste(
              "Peak Exon Type is a positional composition view. It shows whether",
              "saved peaks are more often associated with first, middle, or last",
              "exon overlaps."
            ),
            paste(
              "This page is usually the quickest exon-level summary to read because",
              "it emphasizes composition rather than a continuous distribution."
            )
          )
        ),
        list(
          title = "Peak Exon Num",
          paragraphs = c(
            paste(
              "Peak Exon Num summarizes how many exons belong to the transcripts",
              "associated with the saved peaks. It helps distinguish simpler from",
              "more exon-rich transcript structures."
            ),
            paste(
              "Read it together with Peak Exon Size when you want both length and",
              "structural-complexity context."
            )
          )
        )
      )
    ),
    list(
      id = "site",
      label = "Site",
      summary = paste(
        "Profile local signal around transcription, translation, and splice-site",
        "boundaries."
      ),
      blocks = list(
        list(
          title = "Transcription",
          paragraphs = c(
            paste(
              "Transcription profiles density around transcription start and end",
              "neighborhoods. In RNAmeta Online it provides both density curves",
              "and heatmap panels."
            ),
            paste(
              "When multiple files are saved, the heatmap opens on the first saved",
              "sample by default and can be switched to another saved sample."
            )
          )
        ),
        list(
          title = "Translation",
          paragraphs = c(
            paste(
              "Translation focuses on neighborhoods around translation start and",
              "end boundaries. Like Transcription, it combines density curves with",
              "a heatmap view tied to the current saved session."
            ),
            paste(
              "Use it when your interpretation depends more on coding-boundary",
              "behavior than on whole-transcript structure."
            )
          )
        ),
        list(
          title = "Splicesite",
          paragraphs = c(
            paste(
              "Splicesite measures density around 5' and 3' splice-site regions.",
              "Unlike Transcription and Translation, this view is focused on",
              "density profiles rather than a density-plus-heatmap display."
            ),
            paste(
              "This is the most local and junction-oriented module in the Site",
              "group."
            )
          )
        )
      )
    ),
    list(
      id = "motif",
      label = "Motif",
      summary = paste(
        "Run motif analysis in the background, monitor job progress, and",
        "reopen completed results by job id."
      ),
      blocks = list(
        list(
          title = "What the Motif page does",
          paragraphs = c(
            paste(
              "Motif is the web-specific analysis workflow for motif discovery in",
              "the current RNAmeta Online release. It runs HOMER, MEME, or STREME in the",
              "background and lets you return later with the saved job id to review",
              "the resulting motif patterns."
            ),
            paste(
              "The page reuses the BED files saved in Upload / Run. Choose the",
              "saved sample, motif engine, peak subset, motif width, and maximum",
              "motif count before submitting the job."
            )
          ),
          lists = list(
            list(
              label = "Motif workflow checkpoints",
              items = c(
                "Save at least one BED file in Upload / Run before opening Motif.",
                "Watch the queue status card to see how many jobs are active or queued.",
                "Keep the job id if you want to reopen the same result later."
              )
            )
          )
        ),
        list(
          title = "Running and reopening jobs",
          paragraphs = c(
            paste(
              "Submitted jobs continue in the background even if you leave the page.",
              "Use the Motif Job ID field and Query Job to reopen a queued, running,",
              "failed, or completed job."
            ),
            paste(
              "Completed jobs unlock the export panel and render motif output in the",
              "main workspace. If a job completes but RNAmeta Online cannot parse a usable",
              "motif result, try a larger motif width or a different BED file."
            )
          )
        )
      )
    ),
    list(
      id = "export-cache",
      label = "Export & Cache",
      summary = paste(
        "Understand when export panels appear, what figure and table exports",
        "contain, and how runtime caching affects repeated analyses."
      ),
      blocks = list(
        list(
          title = "Figure and table export",
          paragraphs = c(
            paste(
              "Most analysis pages show an export panel only after a successful run.",
              "Figure export is intended for chart capture, while data export writes",
              "the machine-readable values behind the current view."
            ),
            paste(
              "Depending on the module, CSV or TXT export may contain normalized",
              "density tables, feature counts, overlap summaries, or combined",
              "density-plus-heatmap tables."
            )
          ),
          lists = list(
            list(
              label = "Export guidance",
              items = c(
                "Use figure export when the visual chart is the deliverable.",
                "Use CSV or TXT when you need the values behind the current chart state.",
                "If a page lets you hide features or switch samples, export reflects the active view."
              )
            )
          )
        ),
        list(
          title = "What runtime cache means in RNAmeta Online",
          paragraphs = c(
            paste(
              "RNAmeta Online keeps runtime caches for repeated work inside the active",
              "analysis session. This is why rerunning the same inputs can be faster",
              "than the first pass."
            ),
            paste(
              "Treat caching as a speed feature rather than a permanent result store.",
              "For shareable or long-lived outputs, rely on exported figures,",
              "tables, or saved motif job ids instead of cache state."
            )
          )
        )
      )
    ),
    list(
      id = "faq",
      label = "FAQ / Troubleshooting",
      summary = paste(
        "Resolve the most common web workflow problems before assuming the",
        "analysis logic is wrong."
      ),
      blocks = list(
        list(
          title = "Why is the Run button disabled",
          paragraphs = c(
            paste(
              "Most modules require at least one saved BED file in Upload / Run.",
              "Gene Matrix needs 2 to 5 saved files. If the page shows a hint about",
              "saving BED files first, go back to Upload / Run, confirm species,",
              "and click Save."
            )
          )
        ),
        list(
          title = "Why did the result stop matching my new inputs",
          paragraphs = c(
            paste(
              "Changing the saved Upload / Run context means the downstream result",
              "you are looking at may no longer match the current files or species.",
              "Rerun the affected modules after every important input change."
            )
          )
        ),
        list(
          title = "Why does exported data not match hidden features",
          paragraphs = c(
            paste(
              "Export follows the active rendered state. For example, Peak",
              "Distribution exports the feature set that is currently visible, not",
              "every feature that exists behind the page."
            )
          )
        ),
        list(
          title = "Why can I not find a motif result immediately",
          paragraphs = c(
            paste(
              "Motif jobs are queued and run in the background. Check the queue",
              "status card, keep the job id, and use Query Job to reopen the result",
              "later if the analysis has not finished yet."
            ),
            paste(
              "If a completed job cannot be parsed into a supported motif result,",
              "try a larger motif width or a BED file with wider intervals."
            )
          )
        )
      )
    )
  )
}
