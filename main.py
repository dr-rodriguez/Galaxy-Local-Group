from astropy.table import Table

def main():
    t = Table.read("NearbyGalaxies_Jan2021_PUBLIC.fits")
    meta = t.meta
    
    output_lines = ["# Galaxy Table Columns\n", "| Column | Data Type | Description |\n", "|--------|-----------|-------------|\n"]
    
    for i, column in enumerate(t.colnames):
        col_id = f"TCOMM{i+1}"
        extra_metadata = meta.get(col_id, "")
        output_lines.append(f"| {column} | {t[column].dtype} | {extra_metadata} |\n")
    
    with open("galaxy_columns.md", "w") as f:
        f.writelines(output_lines)


if __name__ == "__main__":
    main()
