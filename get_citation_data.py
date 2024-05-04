from pypdb import get_info
# print(list(all_info.keys()))

def get_citation_data(pdbid):
    all_info = get_info(pdbid)
    
    if 'citation' not in all_info:
        return "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"
    citation = all_info['citation'][0]
    
    authors = "N/A"
    if 'rcsb_authors' in citation:
        authors = citation['rcsb_authors']
    
    citation_title = "N/A"
    if 'title' in citation:
        citation_title = citation['title']

    doi = "?"
    if 'pdbx_database_id_doi' in citation:
        doi = citation['pdbx_database_id_doi']

    year = "N/A"
    if 'year' in citation:
        year = citation['year']
    
    pubmed = "?"
    if 'pdbx_database_id_pub_med' in citation:
        pubmed = citation['pdbx_database_id_pub_med']

    method = "N/A"
    if 'exptl' in all_info:
        # print("YUP!")
        if 'method' in all_info['exptl'][0]:
            method = all_info['exptl'][0]['method']

    keywords = "N/A"
    if 'struct_keywords' in all_info:
        if 'text' in all_info['struct_keywords']:
            keywords = all_info['struct_keywords']['text']

    release_date = "?"
    if 'rcsb_accession_info' in all_info:
        if 'initial_release_date' in all_info['rcsb_accession_info']:
            release_date = all_info['rcsb_accession_info']['initial_release_date'][0:10]
    
    title = "N/A"
    if 'struct' in all_info:
        if 'title' in all_info['struct']:
            title = all_info['struct']['title']

    return citation_title, year, authors, doi, pubmed, method, keywords, release_date, title

if __name__ == "__main__":
    get_citation_data('1jgg')