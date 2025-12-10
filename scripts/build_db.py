"""
Build the database from the JSON files.
"""

from astrodb_utils import build_db_from_json
from astrodb_utils.loaders import DatabaseSettings

def main():
    # Instantiate database from settings file
    db_settings = DatabaseSettings(settings_file="database.toml")
    _ = build_db_from_json(settings_file=db_settings.settings_file)

if __name__ == "__main__":
    main()
