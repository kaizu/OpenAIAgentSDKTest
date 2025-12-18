# -*- coding: utf-8 -*-
"""Entrez agent tools."""

import logging
import os
from typing import TypedDict
from agents import Agent, function_tool
from utils import my_custom_error_function
from Bio import Entrez, Medline

logger = logging.getLogger(__name__)

def _ensure_entrez_email() -> str:
    if not Entrez.email:
        email = os.getenv("ENTREZ_EMAIL")
        if not email:
            raise RuntimeError("ENTREZ_EMAIL が環境変数に設定されていません。")
        Entrez.email = email
    return Entrez.email

@function_tool(failure_error_function=my_custom_error_function)
async def search_pubmed(term: str, retmax: int = 10) -> list[str]:
    """PubMedで単語の検索を行い、文献のIDのリストを返します。
    
    Args:
        term: 検索したい単語。
        retmax: 返す文献の数の最大値。0より大きく、100以下の整数で与えること。
    """
    logger.info(f"search_pubmed called with term={term}, retmax={retmax}")
    assert 0 < retmax <= 100
    _ensure_entrez_email()
    handle = Entrez.esearch(db="pubmed", term=term, retmax=f"{retmax:d}")
    record = Entrez.read(handle)
    return record["IdList"][: 10]

class PubMedDocumentInfo(TypedDict):
    id: str
    title: str
    author: list[str]
    pubdate: str
    source: str

@function_tool(failure_error_function=my_custom_error_function)
async def fetch_pubmed_document_summary(primary_ids: list[str]) -> list[PubMedDocumentInfo]:
    """PubMedの文献のIDのリストから、それぞれの情報をリストにして返します。
    この文献情報はid, title, author, pubdate, sourceからなり、
    それぞれ文献ID、タイトル、著者名を要素とするリスト、発表年月日、雑誌名です。
    
    Args:
        primary_ids: 情報を取得したいIDのリスト
    """
    logger.info(f"fetch_pubmed_document_summary called for ids={primary_ids}")
    assert len(primary_ids) <= 100
    _ensure_entrez_email()
    handle = Entrez.esummary(db="pubmed", id=primary_ids)
    record = Entrez.read(handle)
    result = [
        PubMedDocumentInfo(id=r['Id'], title=r["Title"], author=r['AuthorList'], pubdate=r['PubDate'], source=r['Source'])
        for r in record]
    return result

@function_tool(failure_error_function=my_custom_error_function)
async def fetch_pubmed_document_abstract(primary_ids: list[str]) -> list[str]:
    """PubMedの文献のIDのリストから、それぞれの要旨をリストにして返します。
    
    Args:
        primary_ids: 要旨を取得したい文献のIDのリスト
    """
    logger.info(f"fetch_pubmed_document_abstract called for ids={primary_ids}")
    assert len(primary_ids) <= 100
    _ensure_entrez_email()
    handle = Entrez.efetch(db="pubmed", id=primary_ids, rettype="medline", retmode="text")
    record = Medline.parse(handle)
    result = [r.get('AB', '') for r in record]
    return result

def create_entrez_agent(model: str) -> Agent:
    logger.info("create_entrez_agent called")
    return Agent(
            name="Entrez assistant agent",
            model=model,
            tools=[
                search_pubmed,
                fetch_pubmed_document_summary,
                fetch_pubmed_document_abstract,
                ],
            instructions=(
                "あなたはPubMed上の文献の検索を支援するアシスタントです。"
                "ツールを使って文献を検索し、得られた結果を要約して返してください。"
            ),
        )
