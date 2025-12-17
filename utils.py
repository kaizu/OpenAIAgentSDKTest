# -*- coding: utf-8 -*-
"""Shared utilities."""

import logging
from typing import Any
from agents import RunContextWrapper

logger = logging.getLogger(__name__)


def my_custom_error_function(context: RunContextWrapper[Any], error: Exception) -> str:
    """Provide a user-friendly error message for tool failures."""
    logger.info(f"my_custom_error_function called with error: {error}")
    return f"ツール実行でエラーが発生しました: {error}"
