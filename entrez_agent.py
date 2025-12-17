import logging
from typing import TypedDict
from agents import Agent, function_tool
from utils import my_custom_error_function
from Bio import Entrez, Medline

logger = logging.getLogger(__name__)

Entrez.email = "kwaizu@gmail.com"     # Always tell NCBI who you are

@function_tool(failure_error_function=my_custom_error_function)
async def search_pubmed(term: str) -> list[str]:
    """PubMedで単語の検索を行い、文献のIDのリストを返します。
    
    Args:
        term: 検索したい単語。
    """
    logger.info(f"search_pubmed called with term={term}")
    handle = Entrez.esearch(db="pubmed", term=term, retmax="10")
    record = Entrez.read(handle)
    return record["IdList"][: 10]

class PubMedDocumentInfo(TypedDict):
    id: str
    title: str
    author: list[str]

@function_tool(failure_error_function=my_custom_error_function)
async def fetch_pubmed_document_summary(primary_ids: list[str]) -> list[PubMedDocumentInfo]:
    """PubMedの文献のIDのリストから、それぞれの情報をリストにして返します。
    この文献情報はid, title, authorからなり、それぞれ文献ID、タイトル、著者名を要素とするリストです。
    
    Args:
        primary_ids: 情報を取得したいIDのリスト
    """
    handle = Entrez.esummary(db="pubmed", id=primary_ids)
    record = Entrez.read(handle)
    result = [PubMedDocumentInfo(id=r['Id'], title=r["Title"], author=r['AuthorList']) for r in record]
    return result

@function_tool(failure_error_function=my_custom_error_function)
async def fetch_pubmed_document_abstract(primary_ids: list[str]) -> list[str]:
    """PubMedの文献のIDのリストから、それぞれの要旨をリストにして返します。
    
    Args:
        primary_ids: 要旨を取得したい文献のIDのリスト
    """
    handle = Entrez.efetch(db="pubmed", id=primary_ids, rettype="medline",retmode="text")
    record = Medline.parse(handle)
    result = [r.get('AB', '') for r in record]
    return result

def create_entrez_agent() -> Agent:
    logger.info("create_entrez_agent called")
    return Agent(
            name="EntrezAssistant",
            # model="gpt-4o-mini",
            model="gpt-5.2",
            tools=[
                search_pubmed,
                fetch_pubmed_document_summary,
                fetch_pubmed_document_abstract,
                ],
            instructions=(
                "あなたはNCBI Entrez検索を支援するアシスタントです。"
                "ツールを使って検索し、得られた結果を要約して返してください。"
            ),
        )
