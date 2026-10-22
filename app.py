"""
========================================================================================
PubMed & Gene Explorer - Bioinformatics Research Web Application
========================================================================================
A beginner-friendly web application designed for BSc Biotechnology students to explore
scientific publications, analyze gene-disease literature, visualize publication trends,
and generate evidence-based research summary reports.

Author: Senior Python Developer & Bioinformatics Specialist
Tech Stack: Python 3.10+, Streamlit, Biopython (NCBI Entrez), Pandas, Plotly
========================================================================================
"""

import re
import datetime
from collections import Counter
from typing import List, Dict, Tuple, Optional, Any

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio import Entrez
import urllib.error


# ======================================================================================
# 1. PAGE CONFIGURATION & CUSTOM STYLING
# ======================================================================================

st.set_page_config(
    page_title="PubMed & Gene Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for scientific/biotechnology aesthetic and polished card designs
st.markdown(
    """
    <style>
    /* Main container styling */
    .main {
        background-color: #f8fafc;
    }
    
    /* Sidebar container styling & width constraints */
    [data-testid="stSidebar"], section[data-testid="stSidebar"] {
        min-width: 320px !important;
        max-width: 380px !important;
        width: 340px !important;
    }
    
    /* Prevent sidebar horizontal scrolling */
    [data-testid="stSidebar"] > div:first-child,
    [data-testid="stSidebarContent"],
    div[data-testid="stSidebarUserContent"] {
        overflow-x: hidden !important;
        box-sizing: border-box !important;
    }
    
    /* Preset grid buttons: single line, neat alignment, no text wrapping */
    [data-testid="stSidebar"] .stButton > button {
        white-space: nowrap !important;
        text-overflow: ellipsis !important;
        padding: 0.35rem 0.3rem !important;
        font-size: 0.85rem !important;
        font-weight: 600 !important;
        min-height: 2.2rem !important;
        height: 2.2rem !important;
        width: 100% !important;
        display: flex !important;
        align-items: center !important;
        justify-content: center !important;
    }
    [data-testid="stSidebar"] .stButton > button p {
        white-space: nowrap !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;
        margin: 0 !important;
        font-size: 0.85rem !important;
    }

    /* Primary search button in sidebar */
    [data-testid="stSidebar"] button[kind="primary"],
    [data-testid="stSidebar"] .stButton > button[data-testid="baseButton-primary"] {
        min-height: 2.6rem !important;
        height: 2.6rem !important;
        font-size: 0.95rem !important;
        font-weight: 700 !important;
    }
    
    /* Column blocks in sidebar: tidy spacing, no overflow */
    [data-testid="stSidebar"] [data-testid="stHorizontalBlock"] {
        gap: 0.45rem !important;
        width: 100% !important;
    }

    /* Uniform alignment for input fields, dropdowns, and sliders */
    [data-testid="stSidebar"] .stTextInput,
    [data-testid="stSidebar"] .stSelectbox,
    [data-testid="stSidebar"] .stSlider {
        width: 100% !important;
        margin-bottom: 0.35rem !important;
    }
    
    /* Header branding */
    .app-title {
        font-family: 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
        color: #0f172a;
        font-weight: 800;
        font-size: 2.2rem;
        margin-bottom: 0.2rem;
        display: flex;
        align-items: center;
        gap: 0.5rem;
    }
    .app-subtitle {
        color: #475569;
        font-size: 1.05rem;
        margin-bottom: 1.5rem;
    }
    
    /* Paper card container */
    .paper-card {
        background-color: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 1.25rem 1.5rem;
        margin-bottom: 1.25rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        transition: transform 0.15s ease, box-shadow 0.15s ease;
    }
    .paper-card:hover {
        border-color: #0284c7;
        box-shadow: 0 4px 12px rgba(2, 132, 199, 0.08);
    }
    
    /* Paper title link */
    .paper-title {
        font-size: 1.15rem;
        font-weight: 700;
        color: #0369a1 !important;
        text-decoration: none;
        line-height: 1.4;
    }
    .paper-title:hover {
        text-decoration: underline;
        color: #0284c7 !important;
    }
    
    /* Metadata pills */
    .badge-pill {
        display: inline-block;
        padding: 0.2rem 0.55rem;
        font-size: 0.78rem;
        font-weight: 600;
        border-radius: 9999px;
        margin-right: 0.4rem;
        margin-top: 0.3rem;
        margin-bottom: 0.5rem;
    }
    .badge-pmid {
        background-color: #e0f2fe;
        color: #0369a1;
        border: 1px solid #bae6fd;
    }
    .badge-year {
        background-color: #f1f5f9;
        color: #334155;
        border: 1px solid #cbd5e1;
    }
    .badge-journal {
        background-color: #fef3c7;
        color: #92400e;
        border: 1px solid #fde68a;
    }
    
    /* Abstract container */
    .abstract-box {
        background-color: #f8fafc;
        border-left: 3px solid #0284c7;
        padding: 0.75rem 1rem;
        border-radius: 4px;
        font-size: 0.92rem;
        color: #334155;
        line-height: 1.6;
        margin-top: 0.5rem;
    }
    
    /* Disclaimer badge */
    .evidence-disclaimer {
        background-color: #eff6ff;
        border-left: 4px solid #3b82f6;
        padding: 0.8rem 1rem;
        border-radius: 6px;
        font-size: 0.88rem;
        color: #1e3a8a;
        margin-bottom: 1.5rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# ======================================================================================
# 2. QUERY VALIDATION & FORMULATION FUNCTIONS
# ======================================================================================

def validate_query(query: str) -> Tuple[bool, str]:
    """
    Validates user input to prevent empty searches or invalid punctuation-only queries.
    
    Parameters:
        query (str): The search string entered by the user.
        
    Returns:
        Tuple[bool, str]: (is_valid, error_or_warning_message)
    """
    if not query:
        return False, "Search query cannot be empty. Please enter a gene or disease name."
    
    stripped = query.strip()
    if len(stripped) == 0:
        return False, "Search query cannot be just whitespace. Please enter a gene or disease name."
    
    # Check if the query contains at least one alphanumeric character
    if not re.search(r"[A-Za-z0-9]", stripped):
        return False, "Search query must contain at least one letter or number (e.g., 'BRCA1', 'Diabetes')."
    
    if len(stripped) < 2:
        return False, "Query is too short. Please enter at least 2 characters for meaningful results."
        
    return True, ""


def build_pubmed_query(term: str, search_mode: str) -> str:
    """
    Constructs an optimized NCBI Entrez search syntax based on the selected mode.
    
    Bioinformatics Note:
    - [Gene]: Restricts search to official gene symbols/names in Entrez.
    - [Title/Abstract]: Searches paper title and abstract fields.
    - [MeSH Terms]: Medical Subject Headings controlled vocabulary thesaurus.
    """
    term_clean = term.strip()
    
    if search_mode == "Gene Focus":
        # Search both the specific gene tag and title/abstract for comprehensive coverage
        return f"({term_clean}[Gene Name]) OR ({term_clean}[Title/Abstract])"
    elif search_mode == "Disease Focus":
        # Search disease in Title/Abstract and MeSH terms
        return f"({term_clean}[Title/Abstract]) OR ({term_clean}[MeSH Terms])"
    elif search_mode == "Exact Phrase":
        return f'"{term_clean}"'
    else:
        # Default: All Fields
        return term_clean


# ======================================================================================
# 3. LIVE PUBMED DATA RETRIEVAL (NCBI ENTREZ API)
# ======================================================================================

@st.cache_data(show_spinner=False, ttl=1800)
def fetch_pubmed_data(
    query: str,
    max_results: int = 25,
    sort_by: str = "pub_date",
    email: str = "student.biotech@example.com",
    api_key: Optional[str] = None,
) -> Tuple[pd.DataFrame, Optional[str]]:
    """
    Fetches real-time scientific articles from PubMed using the NCBI Entrez API via Biopython.
    
    How this works:
    1. Entrez.esearch: Searches PubMed index with query, returns a list of PubMed IDs (PMIDs).
    2. Entrez.efetch: Downloads full XML records for those PMIDs in a single batch.
    3. XML Parsing: Extracts Title, Abstract, Year, Authors, Journal, and builds links.
    
    Parameters:
        query (str): PubMed query string.
        max_results (int): Number of papers to retrieve (10 to 100).
        sort_by (str): 'pub_date' (Most Recent) or 'relevance' (Best Match).
        email (str): Required by NCBI Entrez to notify user in case of issues.
        api_key (str): Optional NCBI API key for higher rate limits (10 req/s vs 3 req/s).
        
    Returns:
        Tuple[pd.DataFrame, Optional[str]]: (DataFrame of articles, error_message)
    """
    # NCBI Entrez requirement: Email must be set
    Entrez.email = email.strip() if email else "student.biotech@example.com"
    Entrez.tool = "PubMedGeneExplorer_BSc"
    if api_key and api_key.strip():
        Entrez.api_key = api_key.strip()
    
    try:
        # Step 1: ESearch - Find PMIDs matching query
        search_handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            sort=sort_by
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        pmid_list = search_results.get("IdList", [])
        
        # If no papers match query
        if not pmid_list:
            return pd.DataFrame(), None
        
        # Step 2: EFetch - Batch download XML records for all retrieved PMIDs
        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=pmid_list,
            rettype="xml",
            retmode="xml"
        )
        articles_xml = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        # Step 3: Parse XML structures into clean structured dictionary records
        records = []
        for article_data in articles_xml.get("PubmedArticle", []):
            medline = article_data.get("MedlineCitation", {})
            article = medline.get("Article", {})
            
            # 3a. PubMed ID (PMID)
            pmid = str(medline.get("PMID", "N/A"))
            
            # 3b. Article Title
            title = article.get("ArticleTitle", "No title available").strip()
            # Clean trailing periods or brackets
            if title.startswith("[") and title.endswith("]"):
                title = title[1:-1]
                
            # 3c. Abstract (handles both simple text and multi-section structured abstracts)
            abstract_text = "No abstract available for this record (citation or title-only entry)."
            abstract_obj = article.get("Abstract", {})
            if "AbstractText" in abstract_obj:
                raw_text = abstract_obj["AbstractText"]
                if isinstance(raw_text, list):
                    # Structured abstracts have sections like BACKGROUND, METHODS, RESULTS, CONCLUSION
                    section_texts = []
                    for section in raw_text:
                        label = getattr(section, "attributes", {}).get("Label", "") if hasattr(section, "attributes") else ""
                        text_val = str(section).strip()
                        if label and label.lower() != "unlabelled":
                            section_texts.append(f"**{label}**: {text_val}")
                        else:
                            section_texts.append(text_val)
                    abstract_text = "\n\n".join(section_texts)
                else:
                    abstract_text = str(raw_text).strip()
            
            # 3d. Authors list
            author_names = []
            for author in article.get("AuthorList", []):
                last_name = author.get("LastName", "")
                fore_name = author.get("ForeName", "") or author.get("Initials", "")
                if last_name:
                    author_names.append(f"{last_name} {fore_name}".strip())
                elif "CollectiveName" in author:
                    author_names.append(author["CollectiveName"].strip())
            
            if author_names:
                authors_str = ", ".join(author_names)
                first_author = author_names[0]
            else:
                authors_str = "Consortium / Authors not specified"
                first_author = "Unknown"
            
            # 3e. Journal name
            journal_info = article.get("Journal", {})
            journal_title = journal_info.get("Title", "Unknown Journal")
            journal_iso = journal_info.get("ISOAbbreviation", journal_title)
            
            # 3f. Publication Year (Safely extract across multiple possible XML date fields)
            pub_date = journal_info.get("JournalIssue", {}).get("PubDate", {})
            year_val = pub_date.get("Year")
            
            # Fallback 1: MedlineDate format e.g., '2023 Nov-Dec'
            if not year_val and "MedlineDate" in pub_date:
                match = re.search(r"\b(19\d\d|20\d\d)\b", str(pub_date["MedlineDate"]))
                if match:
                    year_val = match.group(1)
            
            # Fallback 2: ArticleDate element
            if not year_val:
                for art_date in article.get("ArticleDate", []):
                    if "Year" in art_date:
                        year_val = art_date["Year"]
                        break
            
            # Fallback 3: MedlineCitation DateCompleted or DateRevised
            if not year_val:
                date_completed = medline.get("DateCompleted", {})
                year_val = date_completed.get("Year")
                
            try:
                year_int = int(year_val) if year_val else None
            except (ValueError, TypeError):
                year_int = None
                
            # 3g. Permanent URL to original PubMed page
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            
            # 3h. DOI identifier if available
            doi = "N/A"
            for item in article_data.get("PubmedData", {}).get("ArticleIdList", []):
                if getattr(item, "attributes", {}).get("IdType") == "doi":
                    doi = str(item)
                    break

            records.append({
                "pmid": pmid,
                "title": title,
                "abstract": abstract_text,
                "has_abstract": abstract_text != "No abstract available for this record (citation or title-only entry).",
                "authors": authors_str,
                "authors_list": author_names,
                "first_author": first_author,
                "journal": journal_title,
                "journal_iso": journal_iso,
                "year": year_int,
                "doi": doi,
                "url": pubmed_url
            })
            
        df = pd.DataFrame(records)
        return df, None
        
    except urllib.error.HTTPError as e:
        if e.code == 429:
            return pd.DataFrame(), "NCBI rate limit exceeded (HTTP 429). Please wait a few seconds and try again, or provide an NCBI API key in the sidebar."
        return pd.DataFrame(), f"NCBI Entrez HTTP Error: {e.code} - {e.reason}. Please check your query syntax."
    except urllib.error.URLError as e:
        return pd.DataFrame(), f"Network connection error: Unable to connect to NCBI PubMed servers ({e.reason}). Check your internet connection."
    except Exception as e:
        return pd.DataFrame(), f"Error retrieving data from NCBI Entrez API: {str(e)}"


# ======================================================================================
# 4. NATURAL LANGUAGE & BIOMEDICAL CONCEPT EXTRACTION
# ======================================================================================

# Common biomedical and general English stopwords to remove during keyword frequency analysis
BIOMEDICAL_STOPWORDS = {
    # Standard English stopwords
    "a", "about", "above", "after", "again", "against", "all", "am", "an", "and", "any", "are", 
    "aren't", "as", "at", "be", "because", "been", "before", "being", "below", "between", "both", 
    "but", "by", "can't", "cannot", "could", "couldn't", "did", "didn't", "do", "does", "doesn't", 
    "doing", "don't", "down", "during", "each", "few", "for", "from", "further", "had", "hadn't", 
    "has", "hasn't", "have", "haven't", "having", "he", "he'd", "he'll", "he's", "her", "here", 
    "here's", "hers", "herself", "him", "himself", "his", "how", "how's", "i", "i'd", "i'll", 
    "i'm", "i've", "if", "in", "into", "is", "isn't", "it", "it's", "its", "itself", "let's", 
    "me", "more", "most", "mustn't", "my", "myself", "no", "nor", "not", "of", "off", "on", 
    "once", "only", "or", "other", "ought", "our", "ours", "ourselves", "out", "over", "own", 
    "same", "shan't", "she", "she'd", "she'll", "she's", "should", "shouldn't", "so", "some", 
    "such", "than", "that", "that's", "the", "their", "theirs", "them", "themselves", "then", 
    "there", "there's", "these", "they", "they'd", "they'll", "they're", "they've", "this", 
    "those", "through", "to", "too", "under", "until", "up", "very", "was", "wasn't", "we", 
    "we'd", "we'll", "we're", "we've", "were", "weren't", "what", "what's", "when", "when's", 
    "where", "where's", "which", "while", "who", "who's", "whom", "why", "why's", "with", 
    "won't", "would", "wouldn't", "you", "you'd", "you'll", "you're", "you've", "your", "yours", 
    "yourself", "yourselves",
    # Academic/biomedical filler words
    "study", "studies", "result", "results", "method", "methods", "conclusion", "conclusions",
    "background", "objective", "aim", "aims", "using", "used", "patient", "patients",
    "group", "groups", "analysis", "analyses", "data", "significant", "significantly",
    "compared", "showed", "associated", "association", "increased", "decreased", "higher",
    "lower", "observed", "role", "roles", "effect", "effects", "human", "humans", "also",
    "may", "can", "well", "including", "present", "performed", "found", "based", "levels",
    "level", "total", "among", "case", "cases", "two", "one", "three", "first", "new"
}

def extract_top_biological_terms(df: pd.DataFrame, top_n: int = 15) -> List[Tuple[str, int]]:
    """
    Extracts the most frequent biological terms and keywords from paper titles and abstracts.
    Demonstrates basic bioinformatics NLP (tokenization, stopword removal, frequency counting).
    """
    if df.empty:
        return []
        
    corpus = " ".join(df["title"].fillna("") + " " + df["abstract"].fillna(""))
    # Tokenize: words with 3 or more letters/hyphens
    words = re.findall(r"\b[A-Za-z][A-Za-z0-9_-]{2,}\b", corpus.lower())
    
    # Filter stopwords and numeric-only strings
    filtered_words = [
        w for w in words 
        if w not in BIOMEDICAL_STOPWORDS and not w.isdigit() and len(w) > 2
    ]
    
    counter = Counter(filtered_words)
    return counter.most_common(top_n)


# ======================================================================================
# 5. INTERACTIVE PLOTLY VISUALIZATIONS (DARK CHARCOAL THEME)
# ======================================================================================

# Shared dark charcoal theme constants for high contrast and visual clarity
DARK_PAPER_BG = "#1e222d"      # Dark charcoal container background
DARK_PLOT_BG = "#161a23"       # Deeper charcoal for inner plotting grid
CHART_GRID_COLOR = "#334155"   # Subtle but clearly visible slate gridlines
CHART_LINE_COLOR = "#475569"   # Visible zero line and axis baseline
TEXT_TITLE_COLOR = "#ffffff"   # Bright white for chart titles
TEXT_AXIS_COLOR = "#f1f5f9"    # High-contrast bright text for axis titles
TEXT_TICK_COLOR = "#cbd5e1"    # Readable light gray for tick labels
TEXT_LABEL_COLOR = "#ffffff"   # Bright white for data bar numbers


def create_publication_trend_chart(df: pd.DataFrame) -> go.Figure:
    """
    Creates an interactive Plotly bar chart displaying the publication count per year.
    Features:
    - High-contrast dark charcoal theme with subtle gridlines
    - Bold, vibrant electric cyan/teal gradient with crisp bar outlines
    - Large, high-contrast readable typography for titles, axes, and outside data labels
    - Responsive layout with unclipped labels and interactive hover tooltips
    """
    valid_years_df = df.dropna(subset=["year"]).copy()
    valid_years_df["year"] = valid_years_df["year"].astype(int)
    
    if valid_years_df.empty:
        # Return empty placeholder figure with dark charcoal theme
        fig = go.Figure()
        fig.add_annotation(
            text="No publication year data available to plot.",
            showarrow=False,
            font=dict(size=15, color="#94a3b8")
        )
        fig.update_layout(
            paper_bgcolor=DARK_PAPER_BG,
            plot_bgcolor=DARK_PLOT_BG,
            font=dict(family="Segoe UI, Roboto, sans-serif", color=TEXT_AXIS_COLOR),
            margin=dict(l=40, r=40, t=60, b=40),
            autosize=True
        )
        return fig
        
    year_counts = valid_years_df["year"].value_counts().reset_index()
    year_counts.columns = ["Year", "Publication Count"]
    year_counts = year_counts.sort_values("Year")
    
    # Bold, vibrant cyan/teal gradient for publication timeline
    fig = px.bar(
        year_counts,
        x="Year",
        y="Publication Count",
        text="Publication Count",
        title="Distribution of Retrieved Publications by Year",
        labels={"Year": "Publication Year", "Publication Count": "Number of Papers"},
        color="Publication Count",
        color_continuous_scale=["#0284c7", "#06b6d4", "#14b8a6", "#2dd4bf"],
    )
    
    fig.update_traces(
        textposition="outside",
        cliponaxis=False,
        textfont=dict(size=13, color=TEXT_LABEL_COLOR, family="Segoe UI, Roboto, sans-serif"),
        hovertemplate="<b>Year: %{x}</b><br>Publications: %{y}<extra></extra>",
        marker=dict(line=dict(width=1.5, color="#38bdf8"))
    )
    
    fig.update_layout(
        title=dict(
            text="<b>Distribution of Retrieved Publications by Year</b>",
            font=dict(size=18, color=TEXT_TITLE_COLOR, family="Segoe UI, Roboto, sans-serif"),
            x=0.02,
            y=0.96
        ),
        paper_bgcolor=DARK_PAPER_BG,
        plot_bgcolor=DARK_PLOT_BG,
        font=dict(family="Segoe UI, Roboto, sans-serif", color=TEXT_AXIS_COLOR),
        margin=dict(l=50, r=40, t=70, b=50),
        coloraxis_showscale=False,
        xaxis=dict(
            title=dict(text="Publication Year", font=dict(size=14, color=TEXT_AXIS_COLOR)),
            tickfont=dict(size=12, color=TEXT_TICK_COLOR),
            gridcolor=CHART_GRID_COLOR,
            zerolinecolor=CHART_LINE_COLOR,
            showline=True,
            linecolor=CHART_LINE_COLOR,
            tickmode="linear",
            dtick=1
        ),
        yaxis=dict(
            title=dict(text="Number of Papers", font=dict(size=14, color=TEXT_AXIS_COLOR)),
            tickfont=dict(size=12, color=TEXT_TICK_COLOR),
            gridcolor=CHART_GRID_COLOR,
            zerolinecolor=CHART_LINE_COLOR,
            showline=True,
            linecolor=CHART_LINE_COLOR,
        ),
        autosize=True
    )
    
    return fig


def create_top_journals_chart(df: pd.DataFrame, top_n: int = 8) -> go.Figure:
    """
    Creates an interactive horizontal bar chart displaying the top publishing journals.
    Features:
    - High-contrast dark charcoal theme with subtle gridlines
    - Bold, vibrant electric royal blue/indigo gradient
    - High-contrast readable typography and unclipped outside data labels
    """
    if df.empty or "journal_iso" not in df.columns:
        fig = go.Figure()
        fig.add_annotation(
            text="No journal data available to plot.",
            showarrow=False,
            font=dict(size=15, color="#94a3b8")
        )
        fig.update_layout(
            paper_bgcolor=DARK_PAPER_BG,
            plot_bgcolor=DARK_PLOT_BG,
            font=dict(family="Segoe UI, Roboto, sans-serif", color=TEXT_AXIS_COLOR),
            margin=dict(l=40, r=40, t=60, b=40),
            autosize=True
        )
        return fig
        
    journal_counts = df["journal_iso"].value_counts().head(top_n).reset_index()
    journal_counts.columns = ["Journal", "Count"]
    journal_counts = journal_counts.sort_values("Count", ascending=True)
    
    # Bold, vibrant royal blue to indigo gradient
    fig = px.bar(
        journal_counts,
        x="Count",
        y="Journal",
        orientation="h",
        text="Count",
        title=f"Top {len(journal_counts)} Publishing Journals",
        labels={"Count": "Articles", "Journal": "Journal Name"},
        color="Count",
        color_continuous_scale=["#2563eb", "#3b82f6", "#60a5fa", "#93c5fd"],
    )
    
    fig.update_traces(
        textposition="outside",
        cliponaxis=False,
        textfont=dict(size=13, color=TEXT_LABEL_COLOR, family="Segoe UI, Roboto, sans-serif"),
        hovertemplate="<b>%{y}</b><br>Articles: %{x}<extra></extra>",
        marker=dict(line=dict(width=1.5, color="#60a5fa"))
    )
    
    fig.update_layout(
        title=dict(
            text=f"<b>Top {len(journal_counts)} Publishing Journals</b>",
            font=dict(size=18, color=TEXT_TITLE_COLOR, family="Segoe UI, Roboto, sans-serif"),
            x=0.02,
            y=0.96
        ),
        paper_bgcolor=DARK_PAPER_BG,
        plot_bgcolor=DARK_PLOT_BG,
        font=dict(family="Segoe UI, Roboto, sans-serif", color=TEXT_AXIS_COLOR),
        margin=dict(l=50, r=40, t=70, b=50),
        coloraxis_showscale=False,
        xaxis=dict(
            title=dict(text="Number of Articles", font=dict(size=14, color=TEXT_AXIS_COLOR)),
            tickfont=dict(size=12, color=TEXT_TICK_COLOR),
            gridcolor=CHART_GRID_COLOR,
            zerolinecolor=CHART_LINE_COLOR,
            showline=True,
            linecolor=CHART_LINE_COLOR,
        ),
        yaxis=dict(
            title=dict(text="Journal Name", font=dict(size=14, color=TEXT_AXIS_COLOR)),
            tickfont=dict(size=12, color=TEXT_TICK_COLOR),
            gridcolor=CHART_GRID_COLOR,
            zerolinecolor=CHART_LINE_COLOR,
            showline=True,
            linecolor=CHART_LINE_COLOR,
        ),
        autosize=True
    )
    
    return fig


def create_keywords_chart(top_terms: List[Tuple[str, int]]) -> go.Figure:
    """
    Creates an interactive bar chart of the most frequent biological keywords extracted.
    Features:
    - High-contrast dark charcoal theme with subtle gridlines
    - Bold, vibrant emerald green gradient
    - High-contrast readable typography and unclipped outside data labels
    """
    if not top_terms:
        fig = go.Figure()
        fig.add_annotation(
            text="No keywords available to plot.",
            showarrow=False,
            font=dict(size=15, color="#94a3b8")
        )
        fig.update_layout(
            paper_bgcolor=DARK_PAPER_BG,
            plot_bgcolor=DARK_PLOT_BG,
            font=dict(family="Segoe UI, Roboto, sans-serif", color=TEXT_AXIS_COLOR),
            margin=dict(l=40, r=40, t=60, b=40),
            autosize=True
        )
        return fig
        
    terms_df = pd.DataFrame(top_terms, columns=["Keyword", "Frequency"])
    terms_df = terms_df.sort_values("Frequency", ascending=True)
    
    # Bold, vibrant emerald green gradient
    fig = px.bar(
        terms_df,
        x="Frequency",
        y="Keyword",
        orientation="h",
        text="Frequency",
        title="Top Recurrent Biological & Clinical Concepts",
        labels={"Frequency": "Occurrence Frequency", "Keyword": "Biological Concept"},
        color="Frequency",
        color_continuous_scale=["#059669", "#10b981", "#34d399", "#6ee7b7"],
    )
    
    fig.update_traces(
        textposition="outside",
        cliponaxis=False,
        textfont=dict(size=13, color=TEXT_LABEL_COLOR, family="Segoe UI, Roboto, sans-serif"),
        hovertemplate="<b>Term: %{y}</b><br>Frequency: %{x}<extra></extra>",
        marker=dict(line=dict(width=1.5, color="#34d399"))
    )
    
    fig.update_layout(
        title=dict(
            text="<b>Top Recurrent Biological & Clinical Concepts</b>",
            font=dict(size=18, color=TEXT_TITLE_COLOR, family="Segoe UI, Roboto, sans-serif"),
            x=0.02,
            y=0.96
        ),
        paper_bgcolor=DARK_PAPER_BG,
        plot_bgcolor=DARK_PLOT_BG,
        font=dict(family="Segoe UI, Roboto, sans-serif", color=TEXT_AXIS_COLOR),
        margin=dict(l=50, r=40, t=70, b=50),
        coloraxis_showscale=False,
        xaxis=dict(
            title=dict(text="Occurrence Frequency", font=dict(size=14, color=TEXT_AXIS_COLOR)),
            tickfont=dict(size=12, color=TEXT_TICK_COLOR),
            gridcolor=CHART_GRID_COLOR,
            zerolinecolor=CHART_LINE_COLOR,
            showline=True,
            linecolor=CHART_LINE_COLOR,
        ),
        yaxis=dict(
            title=dict(text="Biological Concept / Term", font=dict(size=14, color=TEXT_AXIS_COLOR)),
            tickfont=dict(size=12, color=TEXT_TICK_COLOR),
            gridcolor=CHART_GRID_COLOR,
            zerolinecolor=CHART_LINE_COLOR,
            showline=True,
            linecolor=CHART_LINE_COLOR,
        ),
        autosize=True
    )
    
    return fig


# ======================================================================================
# 6. AUTOMATED SUMMARY REPORT GENERATOR (STRICT EVIDENCE DISTINCTION)
# ======================================================================================

def generate_summary_report(
    df: pd.DataFrame,
    query: str,
    search_mode: str,
    top_terms: List[Tuple[str, int]]
) -> str:
    """
    Generates a structured, evidence-grounded summary report in Markdown format.
    
    STRICT SCIENTIFIC INTEGRITY RULES:
    1. Information extracted directly from papers (PMIDs, titles, journals, authors, verbatim findings)
       is clearly separated from synthesized project conclusions.
    2. No citations or findings are invented. Every highlighted claim points to its verified PMID.
    3. Deterministic analysis: Quantifies study design mentions (e.g. clinical, genomic, in vitro).
    """
    if df.empty:
        return "# PubMed & Gene Explorer: Research Report\n\nNo records were retrieved to generate a report."
        
    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    total_papers = len(df)
    valid_years = df["year"].dropna().astype(int)
    min_year = valid_years.min() if not valid_years.empty else "N/A"
    max_year = valid_years.max() if not valid_years.empty else "N/A"
    
    # Top Journals
    top_journals = df["journal_iso"].value_counts().head(3)
    journal_summary_str = ", ".join([f"*{j}* ({c})" for j, c in top_journals.items()])
    
    # Study Design Frequency (deterministic keyword scan in abstracts)
    abstracts_all = " ".join(df["abstract"].fillna("").str.lower())
    study_signals = {
        "Genomic / Sequencing": len(re.findall(r"\b(sequencing|ngs|mutation|variant|genome|dna)\b", abstracts_all)),
        "Expression / Molecular": len(re.findall(r"\b(expression|mrna|protein|pathway|signaling|receptor)\b", abstracts_all)),
        "Clinical / Cohort / Trial": len(re.findall(r"\b(clinical|cohort|trial|patient|prognosis|survival)\b", abstracts_all)),
        "Therapeutic / Drug / Treatment": len(re.findall(r"\b(therapy|treatment|drug|inhibitor|efficacy|response)\b", abstracts_all)),
        "Cellular / In Vitro": len(re.findall(r"\b(cell line|in vitro|culture|assay|transfection)\b", abstracts_all)),
    }
    
    # Build Markdown Document
    md = []
    md.append(f"# 🧬 Research Summary Report: {query}")
    md.append(f"**Generated by PubMed & Gene Explorer (BSc Biotechnology Project)**  ")
    md.append(f"**Report Timestamp:** {now_str} | **Search Mode:** {search_mode}  \n")
    
    md.append("---")
    md.append("### ⚠️ Scientific Evidence & Attribution Standard")
    md.append(
        "> **Notice on Data Integrity:** This report strictly distinguishes **Directly Extracted Evidence** "
        "(sourced verifiably from NCBI PubMed records with PMIDs) from **Synthesized Project Landscape Insights** "
        "(rule-based statistical categorization). No findings, citations, or data points are artificially invented."
    )
    md.append("\n---\n")
    
    # Section 1: Executive Bibliometric Overview
    md.append("## 1. Bibliometric Overview")
    md.append(f"- **Search Query:** `{query}`")
    md.append(f"- **Retrieved Articles Sample:** {total_papers} publications")
    md.append(f"- **Publication Year Range:** {min_year} – {max_year}")
    md.append(f"- **Leading Journals in Sample:** {journal_summary_str}")
    md.append(f"- **Papers with Full Abstracts:** {df['has_abstract'].sum()} of {total_papers}")
    md.append("")
    
    # Section 2: Prominent Biological & Clinical Concepts
    md.append("## 2. Prominent Biological & Clinical Concepts")
    md.append("Top recurrent keywords extracted from paper titles and abstracts (frequency of occurrence):")
    md.append("")
    md.append("| Rank | Biological Concept / Keyword | Occurrence Frequency |")
    md.append("| :--- | :--- | :--- |")
    for idx, (term, count) in enumerate(top_terms[:10], start=1):
        md.append(f"| {idx} | **{term.capitalize()}** | {count} |")
    md.append("")
    
    # Section 3: Directly Extracted Evidence (Paper Highlights with Verifiable PMIDs)
    md.append("## 3. Directly Extracted Evidence from Literature")
    md.append(
        "*The following excerpts and metadata are extracted directly from the official NCBI PubMed database. "
        "Each entry includes the verified PubMed ID (PMID) with an active link to the primary source.*"
    )
    md.append("")
    
    sample_papers = df[df["has_abstract"]].head(5)
    if sample_papers.empty:
        sample_papers = df.head(5)
        
    for _, row in sample_papers.iterrows():
        pmid = row["pmid"]
        title = row["title"]
        authors = row["authors"]
        journal = row["journal_iso"]
        year = row["year"] if pd.notnull(row["year"]) else "N/A"
        url = row["url"]
        
        # Extract a clean, concise excerpt from the abstract
        raw_abst = row["abstract"]
        sentences = [s.strip() for s in re.split(r"(?<=[.!?]) +", raw_abst) if len(s.strip()) > 20]
        
        # Prefer the conclusion or results sentence if structured, else the first 2 salient sentences
        excerpt = " ".join(sentences[:2]) if sentences else raw_abst[:250] + "..."
        if len(excerpt) > 400:
            excerpt = excerpt[:397] + "..."
            
        md.append(f"### • [{title}]({url})")
        md.append(f"- **Verified PMID:** [`{pmid}`]({url}) | **Year:** {year} | **Journal:** *{journal}*")
        md.append(f"- **Authors:** {authors}")
        md.append(f"- **Directly Extracted Evidence Excerpt:**  ")
        md.append(f"  > *\"{excerpt}\"*")
        md.append("")
        
    # Section 4: Synthesized Project Landscape (Deterministic Biological Pattern Analysis)
    md.append("## 4. Synthesized Project Landscape (Pattern Analysis)")
    md.append(
        "*This section represents a deterministic statistical synthesis of methodological signals "
        "identified across the retrieved abstract corpus:*  \n"
    )
    for signal, count in study_signals.items():
        bar = "█" * min(int(count / 2) + 1, 15) if count > 0 else ""
        md.append(f"- **{signal}:** {count} keyword mentions `{bar}`")
        
    md.append("\n**Key Landscape Observations:**")
    top_signal = max(study_signals.items(), key=lambda x: x[1])
    if top_signal[1] > 0:
        md.append(
            f"1. **Primary Focus Area:** The retrieved literature sample is predominantly enriched for **{top_signal[0]}** "
            f"indicators ({top_signal[1]} thematic mentions)."
        )
    md.append(
        f"2. **Chronological Trajectory:** Research on `{query}` shows publication records spanning from "
        f"{min_year} to {max_year}, illustrating continuous scientific interest and evolving methodologies."
    )
    md.append(
        "3. **Biotechnology Project Application:** These findings can serve as foundational background "
        "for literature reviews, seminar presentations, or hypothesis formulation in gene-disease association studies."
    )
    md.append("")
    
    # Section 5: Appendix of All Retrieved Records
    md.append("## 5. Complete Index of Retrieved Publications")
    md.append("| PMID | Year | Journal | Title | Link |")
    md.append("| :--- | :--- | :--- | :--- | :--- |")
    for _, row in df.iterrows():
        y = int(row["year"]) if pd.notnull(row["year"]) else "N/A"
        clean_title = row["title"].replace("|", "/")
        if len(clean_title) > 65:
            clean_title = clean_title[:62] + "..."
        md.append(f"| {row['pmid']} | {y} | {row['journal_iso']} | {clean_title} | [PubMed]({row['url']}) |")
        
    md.append("\n---\n*Report generated using PubMed & Gene Explorer · Powered by NCBI Entrez & Biopython*")
    
    return "\n".join(md)


# ======================================================================================
# 7. USER INTERFACE & APPLICATION FLOW
# ======================================================================================

def main():
    # ----------------------------------------------------------------------------------
    # 7.1 Sidebar: Search Parameters & Controls
    # ----------------------------------------------------------------------------------
    with st.sidebar:
        st.markdown("## 🧬 PubMed & Gene Explorer")
        st.caption("Bioinformatics Project Tool · BSc Biotechnology")
        st.markdown("---")
        
        st.markdown("### 🔍 Search Parameters")
        
        # Preset Quick-Select Examples for Demonstrations & Beginners
        st.markdown("**Quick Example Presets:**")
        preset_cols = st.columns(3, gap="small")
        if preset_cols[0].button("BRCA1", use_container_width=True):
            st.session_state["query_input"] = "BRCA1"
            st.session_state["search_mode_input"] = "Gene Focus"
        if preset_cols[1].button("TP53", use_container_width=True):
            st.session_state["query_input"] = "TP53"
            st.session_state["search_mode_input"] = "Gene Focus"
        if preset_cols[2].button("INS", use_container_width=True):
            st.session_state["query_input"] = "INS"
            st.session_state["search_mode_input"] = "Gene Focus"
            
        preset_cols2 = st.columns(2, gap="small")
        if preset_cols2[0].button("Diabetes", use_container_width=True):
            st.session_state["query_input"] = "Diabetes"
            st.session_state["search_mode_input"] = "Disease Focus"
        if preset_cols2[1].button("Breast Cancer", use_container_width=True):
            st.session_state["query_input"] = "Breast Cancer"
            st.session_state["search_mode_input"] = "Disease Focus"
            
        st.markdown("")
        
        # User Search Input
        user_query = st.text_input(
            "Gene or Disease Name:",
            value=st.session_state.get("query_input", ""),
            placeholder="e.g. BRCA1, TP53, Diabetes, Alzheimer's",
            help="Enter a human or model organism gene symbol, or a medical condition/disease name.",
            key="search_box"
        )
        
        # Search Mode
        search_mode = st.selectbox(
            "Query Search Mode:",
            options=["Gene Focus", "Disease Focus", "All Fields", "Exact Phrase"],
            index=["Gene Focus", "Disease Focus", "All Fields", "Exact Phrase"].index(
                st.session_state.get("search_mode_input", "Gene Focus")
            ),
            help="Gene Focus targets [Gene Name], Disease Focus targets [MeSH Terms] & [Title/Abstract]."
        )
        
        # Maximum papers to fetch
        max_results = st.slider(
            "Max Papers to Retrieve:",
            min_value=10,
            max_value=100,
            value=25,
            step=5,
            help="Higher numbers fetch more publications but take slightly longer."
        )
        
        # Sort Order
        sort_by = st.selectbox(
            "Sort Order:",
            options=["Most Recent", "Relevance (Best Match)"],
            index=0,
            help="Sort by publication date or NCBI relevance algorithm."
        )
        sort_api_param = "pub_date" if sort_by == "Most Recent" else "relevance"
        
        # NCBI Entrez Advanced Settings Expander
        with st.expander("⚙️ NCBI Entrez API Settings"):
            user_email = st.text_input(
                "Entrez Email:",
                value="student.biotech@example.com",
                help="NCBI requires an email to identify the client application."
            )
            api_key = st.text_input(
                "NCBI API Key (Optional):",
                value="",
                type="password",
                help="Optional. Increases NCBI rate limit from 3 to 10 queries per second."
            )
            st.caption("You can obtain a free NCBI API key from your MyNCBI account settings.")
            
        # Search Trigger Button
        search_btn = st.button("🚀 Fetch Live PubMed Data", type="primary", use_container_width=True)
        
        st.markdown("---")
        with st.expander("ℹ️ About This Project"):
            st.markdown(
                """
                **PubMed & Gene Explorer** is an educational bioinformatics tool designed for BSc Biotechnology students.
                
                **Key Technologies:**
                - **NCBI Entrez E-Utilities:** Real-time access to >36 million biomedical citations.
                - **Biopython:** Scientific library parsing Medline XML structures.
                - **Pandas:** Tabular bioinformatics data structures.
                - **Plotly:** Interactive charts with hover and zoom.
                """
            )

    # ----------------------------------------------------------------------------------
    # 7.2 Main Content Area: Header & Branding
    # ----------------------------------------------------------------------------------
    st.markdown('<div class="app-title">🧬 PubMed & Gene Explorer</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="app-subtitle">Real-time biomedical literature retrieval, bibliometric trend visualization, '
        'and evidence-grounded summary generation for bioinformatics research.</div>',
        unsafe_allow_html=True
    )

    # Manage search state in st.session_state
    if "current_df" not in st.session_state:
        st.session_state["current_df"] = None
    if "current_query" not in st.session_state:
        st.session_state["current_query"] = None
    if "current_mode" not in st.session_state:
        st.session_state["current_mode"] = None

    # Handle Search Button Click
    if search_btn:
        is_valid, validation_msg = validate_query(user_query)
        if not is_valid:
            st.error(f"⚠️ Validation Error: {validation_msg}")
        else:
            with st.spinner(f"Querying NCBI PubMed for '{user_query}' ({search_mode})..."):
                formatted_query = build_pubmed_query(user_query, search_mode)
                df, error_msg = fetch_pubmed_data(
                    query=formatted_query,
                    max_results=max_results,
                    sort_by=sort_api_param,
                    email=user_email,
                    api_key=api_key
                )
                
                if error_msg:
                    st.error(f"❌ {error_msg}")
                elif df.empty:
                    st.warning(
                        f"🔍 No publications found matching query: **'{user_query}'**.\n\n"
                        "**Troubleshooting Suggestions:**\n"
                        "- Check the spelling of the gene symbol or disease.\n"
                        "- Try using official HGNC symbols (e.g. *TP53* instead of *p53*, *ERBB2* instead of *HER2*).\n"
                        "- Switch the Search Mode in the sidebar to **'All Fields'** for broader coverage."
                    )
                else:
                    st.session_state["current_df"] = df
                    st.session_state["current_query"] = user_query
                    st.session_state["current_mode"] = search_mode
                    st.success(f"✅ Successfully retrieved **{len(df)}** publications from PubMed!")

    # ----------------------------------------------------------------------------------
    # 7.3 Results Dashboard or Welcome View
    # ----------------------------------------------------------------------------------
    df = st.session_state.get("current_df")
    active_query = st.session_state.get("current_query")
    active_mode = st.session_state.get("current_mode")

    if df is None:
        # Welcome Screen for Beginners
        st.info("👈 Enter a Gene (e.g., **BRCA1**, **TP53**) or Disease (e.g., **Diabetes**) in the sidebar and click **'Fetch Live PubMed Data'** to begin.")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(
                """
                ### 1. 🔍 Live PubMed Search
                - Direct connection to NCBI Entrez E-Utilities.
                - Validated search syntax with Gene and Disease focus modes.
                - Live abstracts, author lists, and PubMed links.
                """
            )
        with col2:
            st.markdown(
                """
                ### 2. 📊 Publication Trends
                - Interactive Plotly distribution by year.
                - Leading journals & active researcher groups.
                - Recurrent biomedical concept frequency.
                """
            )
        with col3:
            st.markdown(
                """
                ### 3. 📝 Automated Summary Report
                - Strictly grounded in retrieved papers.
                - Direct citations with verified PMIDs.
                - Instant Markdown (`.md`) export for project reports.
                """
            )
        return

    # If results exist, display Key Metrics Bar
    valid_years = df["year"].dropna().astype(int)
    year_range_str = f"{valid_years.min()} – {valid_years.max()}" if not valid_years.empty else "N/A"
    top_journal_name = df["journal_iso"].value_counts().index[0] if not df.empty else "N/A"
    
    metric_cols = st.columns(4)
    metric_cols[0].metric("Retrieved Articles", f"{len(df)} papers")
    metric_cols[1].metric("Publication Span", year_range_str)
    metric_cols[2].metric("Top Journal", str(top_journal_name)[:20])
    metric_cols[3].metric("Abstracts Available", f"{df['has_abstract'].sum()} / {len(df)}")

    st.markdown("---")

    # Extract keywords for analysis tabs
    top_terms = extract_top_biological_terms(df, top_n=15)

    # ----------------------------------------------------------------------------------
    # 7.4 Three Required Tabs
    # ----------------------------------------------------------------------------------
    tab1, tab2, tab3 = st.tabs([
        "📄 Search Results & Abstracts",
        "📈 Publication Trends",
        "📝 Automated Summary Report"
    ])

    # ----------------------------------------------------------------------------------
    # TAB 1: Search Results & Abstracts
    # ----------------------------------------------------------------------------------
    with tab1:
        st.subheader(f"Publications Matching '{active_query}'")
        
        # Filter & Sort controls within retrieved results
        filter_col, sort_col, view_col = st.columns([3, 2, 2])
        filter_text = filter_col.text_input(
            "Filter within results:",
            placeholder="Type keyword, author, or journal...",
            label_visibility="collapsed"
        )
        sort_choice = sort_col.selectbox(
            "Sort by:",
            options=["Year (Newest First)", "Year (Oldest First)", "Title (A-Z)"],
            label_visibility="collapsed"
        )
        show_abstracts_default = view_col.checkbox("Expand All Abstracts", value=False)
        
        # Apply filtering
        display_df = df.copy()
        if filter_text:
            ft = filter_text.lower()
            mask = (
                display_df["title"].str.lower().str.contains(ft, na=False) |
                display_df["abstract"].str.lower().str.contains(ft, na=False) |
                display_df["authors"].str.lower().str.contains(ft, na=False) |
                display_df["journal"].str.lower().str.contains(ft, na=False)
            )
            display_df = display_df[mask]
            st.caption(f"Showing {len(display_df)} of {len(df)} papers matching filter '{filter_text}'")

        # Apply sorting
        if sort_choice == "Year (Newest First)":
            display_df = display_df.sort_values("year", ascending=False)
        elif sort_choice == "Year (Oldest First)":
            display_df = display_df.sort_values("year", ascending=True)
        elif sort_choice == "Title (A-Z)":
            display_df = display_df.sort_values("title", ascending=True)

        if display_df.empty:
            st.info("No papers match the filter text.")
        else:
            # Render individual publication cards
            for _, row in display_df.iterrows():
                pmid = row["pmid"]
                title = row["title"]
                year = int(row["year"]) if pd.notnull(row["year"]) else "Unknown Year"
                journal = row["journal_iso"]
                authors = row["authors"]
                abstract = row["abstract"]
                url = row["url"]
                doi = row["doi"]

                with st.container():
                    st.markdown(
                        f"""
                        <div class="paper-card">
                            <a class="paper-title" href="{url}" target="_blank">🔗 {title}</a>
                            <div style="margin-top: 0.3rem;">
                                <span class="badge-pill badge-pmid">PMID: {pmid}</span>
                                <span class="badge-pill badge-year">Year: {year}</span>
                                <span class="badge-pill badge-journal">{journal}</span>
                            </div>
                            <div style="font-size: 0.88rem; color: #475569; margin-top: 0.2rem;">
                                <b>Authors:</b> {authors}
                            </div>
                        </div>
                        """,
                        unsafe_allow_html=True
                    )
                    
                    # Abstract container
                    with st.expander(f"View Abstract (PMID: {pmid})", expanded=show_abstracts_default):
                        st.markdown(f"**Abstract:**\n\n{abstract}")
                        
                        col_link, col_cite = st.columns([1, 4])
                        col_link.markdown(f"[View on PubMed ↗]({url})")
                        if doi != "N/A":
                            col_cite.caption(f"DOI: https://doi.org/{doi}")

        st.markdown("---")
        # Raw Data CSV Export Button
        csv_data = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="📥 Download Results as CSV",
            data=csv_data,
            file_name=f"pubmed_{active_query.replace(' ', '_')}_results.csv",
            mime="text/csv"
        )

    # ----------------------------------------------------------------------------------
    # TAB 2: Publication Trends
    # ----------------------------------------------------------------------------------
    with tab2:
        st.subheader("Publication & Bibliometric Trends")
        st.caption("Interactive charts powered by Plotly. Hover over data bars for exact counts.")
        
        # 1. Publications by Year (Required Feature)
        trend_fig = create_publication_trend_chart(df)
        st.plotly_chart(trend_fig, use_container_width=True, theme=None)
        
        st.markdown("---")
        
        # 2. Additional Insights: Top Journals & Concept Frequency
        trend_col1, trend_col2 = st.columns(2)
        with trend_col1:
            journals_fig = create_top_journals_chart(df, top_n=8)
            st.plotly_chart(journals_fig, use_container_width=True, theme=None)
            
        with trend_col2:
            keywords_fig = create_keywords_chart(top_terms)
            st.plotly_chart(keywords_fig, use_container_width=True, theme=None)

    # ----------------------------------------------------------------------------------
    # TAB 3: Automated Summary Report
    # ----------------------------------------------------------------------------------
    with tab3:
        st.subheader("Automated Research Summary Report")
        
        # Generate the report text
        report_md = generate_summary_report(df, active_query, active_mode, top_terms)
        
        # Report Export Buttons at the top
        export_col1, export_col2, _ = st.columns([2, 2, 4])
        export_col1.download_button(
            label="📥 Download Report (.md)",
            data=report_md.encode("utf-8"),
            file_name=f"pubmed_summary_{active_query.replace(' ', '_')}.md",
            mime="text/markdown",
            type="primary"
        )
        
        csv_export = df.to_csv(index=False).encode("utf-8")
        export_col2.download_button(
            label="📊 Download Dataset (.csv)",
            data=csv_export,
            file_name=f"pubmed_dataset_{active_query.replace(' ', '_')}.csv",
            mime="text/csv"
        )
        
        st.markdown("---")
        
        # Display the formatted report
        st.markdown(report_md)


if __name__ == "__main__":
    main()
