"""Literature verification of predicted ceRNA triplets (miRNA -> mRNA -> lncRNA).

Checks every (gene, lncRNA) pair of a succeeded mir-network job, one triplet per
bridge miRNA, against PubMed:

  * pubmed.py    -- query builders, rate-limited E-utilities client, SQLite cache
  * rollup.py    -- pure roll-up of per-arm tiers into per-triplet rows
  * review.py    -- pure selection of the hits worth reading
  * symbols.py   -- Ensembl lncRNA transcript/gene id -> HGNC symbol
  * papers.py    -- abstract / full-text fetch and sentence extraction
  * pipeline.py  -- the Celery-side driver that writes the job's files
  * llm_review.py -- standalone python3.11 script: reads the queue, asks the
                     model for a verdict per paper (Anthropic SDK via a gateway)

Everything except llm_review.py is standard-library Python 3.6 and imports
nothing from Flask or Celery, so it is unit-testable in isolation.
"""
