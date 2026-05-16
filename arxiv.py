library(reticulate)
py_run_string("
import arxiv
import pandas as pd

def collect_arxiv_diabetes(total_records):
    client = arxiv.Client()
    count = int(total_records)
    
    search = arxiv.Search(
        query = 'abs:diabetes AND (cat:cs.AI OR cat:cs.LG)',
        max_results = count,
        sort_by = arxiv.SortCriterion.Relevance
    )
    
    results_list = []
    for r in client.results(search):
        results_list.append({
            'PMID': r.get_short_id(),
            'DOI': r.entry_id,
            'Title': r.title,
            'Abstract': r.summary
        })
    
    return pd.DataFrame(results_list)
")

arxiv_df_200 <- py$collect_arxiv_diabetes(200L)
write.csv(arxiv_df_200, "arxiv_diabetes_200.csv", row.names = FALSE)
head(arxiv_df_200)
