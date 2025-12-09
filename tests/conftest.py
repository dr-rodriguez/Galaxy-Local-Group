import os
import pytest
import logging

import astrodb_utils
from astrodb_utils import build_db_from_json
from astrodb_utils.loaders import DatabaseSettings

logger = logging.getLogger(__name__)

# Create a fresh template database for the data and integrity tests
@pytest.fixture(scope="session", autouse=True)
def db():
    logger.info(f"Using version {astrodb_utils.__version__} of astrodb_utils")

    db_settings = DatabaseSettings(settings_file="database.toml")
    db = build_db_from_json(settings_file=db_settings.settings_file)

    # Confirm file was created
    assert os.path.exists(
        db_settings.db_name + ".sqlite"
    ), "Database file 'astrodb-template.sqlite' was not created."

    logger.info(
        f"Loaded {db_settings.db_name} database using build_db_from_json function in conftest.py"
    )

    return db
