# 🧬 PubMed & Gene Explorer

An interactive bioinformatics web application designed for BSc Biotechnology students to search scientific literature, retrieve live PubMed research papers via NCBI Entrez, visualize bibliometric publication trends with Plotly, and generate evidence-based research summary reports.

---

## 🌟 Key Features

1. **Gene and Disease Search**
   - Search by official Gene symbol (e.g., `BRCA1`, `TP53`, `INS`, `EGFR`, `APOE`) or Disease name (e.g., `Diabetes`, `Breast Cancer`, `Alzheimer's`).
   - Query validation: Rejects empty strings, whitespace, or punctuation-only inputs with informative error messages.
   - Presets for quick one-click demo testing.
   - Search modes: **Gene Focus** (`[Gene Name]` & `[Title/Abstract]`), **Disease Focus** (`[MeSH Terms]` & `[Title/Abstract]`), **All Fields**, and **Exact Phrase**.

2. **Live PubMed Data Retrieval (NCBI Entrez)**
   - Connects live to the **NCBI E-Utilities API** via **Biopython**.
   - Extracts Article Title, Abstract (including multi-section structured abstracts), Publication Year, Authors, Journal, DOI, and verified PubMed ID (PMID).
   - Generates clickable links directly to the official PubMed record (`https://pubmed.ncbi.nlm.nih.gov/<PMID>/`).
   - Gracefully handles missing abstracts and edge cases.
   - Caches query results with Streamlit `@st.cache_data` to minimize API latency and respect NCBI rate limits.

3. **Interactive 3-Tab Dashboard**
   - **Tab 1: Search Results & Abstracts**
     - Summary metrics (retrieved papers, publication year span, leading journal).
     - Filter retrieved papers dynamically by keyword, author, or journal without making new API calls.
     - Sort by year (newest/oldest) or title.
     - Expandable cards with full abstracts, badges, and external PubMed links.
     - Export retrieved papers as CSV.
   - **Tab 2: Publication Trends**
     - Interactive **Plotly bar chart** showing the distribution of retrieved publications by year with custom hover metrics.
     - Horizontal bar chart of the top publishing journals.
     - Top recurrent biological and clinical keywords chart (biomedical NLP tokenization).
   - **Tab 3: Automated Summary Report**
     - Generates a concise, structured research summary report.
     - **Strict Evidence Standard**: Clearly separates **Directly Extracted Evidence** (linked to verified PMIDs) from **Synthesized Project Landscape Insights** (rule-based pattern analysis). Never invents findings or fake citations.
     - Downloadable as a formatted **Markdown (`.md`)** file or raw dataset (**CSV**).

---

## 📁 Project Structure

```text
PubMed_Gene_Explorer/
├── app.py              # Main Streamlit web application
├── requirements.txt    # Project dependencies
└── README.md           # Documentation and beginner's guide
```

---

## 🛠️ Installation & Setup (Step-by-Step)

### Prerequisites
- **Python 3.10+** (Tested on Python 3.14 on Windows)
- A terminal (PowerShell, Command Prompt, or VS Code Terminal)
- Internet connection (to query the live NCBI PubMed API)

### 1. Open Terminal in the Project Directory
```powershell
cd c:\Users\LENOVO\Documents\PubMed_Gene_Explorer
```

### 2. (Optional but Recommended) Create a Virtual Environment
```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
```
*(If on Linux/macOS, use `source venv/bin/activate`)*

### 3. Install Dependencies
```powershell
pip install -r requirements.txt
```

The installed libraries are:
- `streamlit`: Web framework for interactive data apps.
- `biopython`: Standard bioinformatics library containing the `Bio.Entrez` module.
- `pandas`: Data manipulation and tabular storage.
- `plotly`: Interactive JavaScript-powered charting library.
- `requests`: HTTP networking library.

### 4. Run the Application
```powershell
streamlit run app.py
```

Your default web browser will automatically open:
```text
http://localhost:8501
```

---
🛠️ Technologies Used

* Python — application logic and data processing.
* Streamlit — interactive web interface.
* NCBI PubMed — biomedical research literature source.
* Biopython Entrez or HTTP requests — PubMed data retrieval, depending on the implementation.
* Plotly — interactive charts and visualizations.
* Pandas — data handling, if used in the application.

📖 How to Use the Project

1. Start the application using the installation instructions above.
2. Enter a gene name, such as BRCA1, TP53, or INS, or a disease name, such as Diabetes or Breast Cancer.
3. Select the available search mode.
4. Choose the maximum number of research papers to retrieve.
5. Run the search and review the results.
6. Open an article’s abstract to read more details.
7. Navigate to Publication Trends to explore the available charts.
8. Open Automated Summary Report to generate and download the report.
9. Download the CSV dataset if you want to examine the retrieved records separately.

⚠️ Limitations

* Results depend on the availability of records in PubMed and the application’s search logic.
* Publication charts reflect the data retrieved by the application and may not represent all publications on a topic.
* Automatically generated summaries should be checked against the original research articles.
* The application is intended for research-literature exploration and is not a substitute for professional medical advice.

📜 Disclaimer

This project is intended for educational and research purposes. Please consult the original publications when interpreting scientific findings.

## 🧬 Understanding the Bioinformatics Code (Beginner's Guide)

If you are presenting this project for your BSc Biotechnology course, here is how the core components work under the hood:

### 1. Connecting to NCBI PubMed with Biopython
NCBI provides a free web service called **E-Utilities** (Entrez Programming Utilities). In `app.py`, we use two primary E-Utilities tools via `Bio.Entrez`:

- **`Entrez.esearch`**: Performs a database query on the PubMed index. It returns matching document IDs (PubMed IDs / PMIDs) rather than full articles:
  ```python
  handle = Entrez.esearch(db="pubmed", term="BRCA1[Gene Name]", retmax=25, sort="pub_date")
  record = Entrez.read(handle)
  pmids = record["IdList"]  # e.g., ['38291024', '38289410', ...]
  ```
- **`Entrez.efetch`**: Downloads the complete XML metadata records for all those PMIDs in a single fast batch request:
  ```python
  fetch_handle = Entrez.efetch(db="pubmed", id=pmids, rettype="xml", retmode="xml")
  articles = Entrez.read(fetch_handle)
  ```

### 2. Parsing PubMed XML
PubMed records are hierarchical XML documents. `app.py` extracts:
- `ArticleTitle`: Title of the paper.
- `AbstractText`: Some papers have structured abstracts with section labels (e.g., `BACKGROUND`, `METHODS`, `RESULTS`, `CONCLUSION`). The code iterates over sections and formats them cleanly with Markdown.
- `AuthorList`: Iterates over authors extracting `LastName` and `ForeName`/`Initials`.
- `PubDate`: Extracts the 4-digit publication year across various formats (such as standard year elements or `MedlineDate` strings like `'2023 Nov-Dec'`).

### 3. Text Processing & Biological Entity Analysis
In `extract_top_biological_terms()`, the application:
1. Concatenates the titles and abstracts into a single text corpus.
2. Tokenizes the text into individual words using regular expressions (`\b[A-Za-z][A-Za-z0-9_-]{2,}\b`).
3. Filters the words against a custom biomedical stopword dictionary to exclude common conversational words (`the`, `with`, `for`) as well as generic academic filler (`study`, `patients`, `results`, `methods`).
4. Uses Python's `collections.Counter` to find the most frequent biological concepts (e.g., `mutation`, `expression`, `signaling`, `pathway`).

### 4. Plotly Interactive Visualizations (Dark Charcoal Theme)
Plotly creates reactive JavaScript charts styled with an accessible, high-contrast dark charcoal theme:
- **Dark Charcoal Theme**: Backgrounds set to `#1e222d` and `#161a23` with visible slate gridlines (`#334155`).
- **High-Contrast Typography**: Bold chart titles (`#ffffff`, 18px), prominent axis titles (`#f1f5f9`, 14px), and tick labels (`#cbd5e1`, 12px) for effortless readability without zooming.
- **Vibrant Color Schemes**:
  - `create_publication_trend_chart()`: Electric cyan/teal gradient displaying chronological research intensity.
  - `create_top_journals_chart()`: Electric royal blue/indigo horizontal bar chart showing leading journals.
  - `create_keywords_chart()`: Vivid emerald green horizontal bar chart displaying recurrent biological terms.
- **Unclipped Value Labels**: Numerical data labels rendered outside bars (`cliponaxis=False`) with interactive hover tooltips.

### 5. Automated Summary Report & Scientific Integrity
In scientific research, separating direct factual quotes from speculative interpretation is essential:
- **Direct Evidence**: Verbatim title, publication year, journal, author list, and actual sentences extracted from the abstract, each tagged with its clickable `[PMID: ...]`.
- **Landscape Synthesis**: Rule-based quantification of study methodologies (genomic vs. clinical vs. in vitro) derived purely from keyword occurrences in the abstracts.
- **Zero Hallucination**: No generative model is permitted to invent PMIDs, authors, or experimental findings.

---

## ❓ Troubleshooting

| Issue | Cause | Solution |
| :--- | :--- | :--- |
| **HTTP 429 Too Many Requests** | NCBI rate limit reached (>3 requests per second) | Wait 5 seconds and retry, or add a free NCBI API Key in the sidebar under "NCBI Entrez API Settings". |
| **Zero publications found** | Typo or obscure synonym in search box | Try official HGNC gene symbols (e.g., *TP53* instead of *p53*), or change Search Mode to "All Fields". |
| **Missing abstracts** | Some older papers or conference notes have no abstract in PubMed | The app gracefully flags these as `"No abstract available for this record"` and allows linking to the full PubMed entry. |
