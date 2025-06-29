"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
import pandas as pd
import re
import warnings

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
        ".fq",
        ".fastq"
    )

    def __init__(
        self, 
        required_columns,
        fastq_dir=None,
        clinical=False,
        single_col="single_end",
        **kwargs):
        """
        Initialize the row checker with the expected required columns.

        Args:
            required_columns (set): A set of column names that are required in the input row.
            single_col (str): The name of the new column that will be inserted and records whether the sample contains single- or paired-end sequencing reads (default "single_end").
            fastq_dir (str): The name of the directory that contains the FASTQ files (default None).
            clinical (bool): Whether the samplesheet is a clinical samplesheet (default False).
            **kwargs: Additional keyword arguments.

        """
        super().__init__(**kwargs)
        self.required_columns = required_columns
        self.fastq_dir = fastq_dir
        self.clinical = clinical
        self._single_col = single_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_row(row)
        self._validate_pair(row)
        self._create_fastq_path(row)
        if self.clinical:
            if row["status"].lower() != "discontinued":
                self._seen.add(tuple(row.values()))
                self.modified.append(row)
        else:
            self._seen.add(tuple(row.values()))
            self.modified.append(row)

        print(row)

        if self.clinical:
            try:
                self._check_barcode_format(row)
                # Barcode format is valid, proceed with further processing
            except ValueError as e:
                print("Error: {}".format(str(e)))
                # Handle the error appropriately, such as logging or terminating the process

    def _validate_row(self, row):
        """Assert that the required columns exist in the input row."""
        if not self.required_columns.issubset(row.keys()):
            missing_columns = self.required_columns - set(row.keys())
            raise AssertionError("Missing required columns: {}".format(", ".join(missing_columns)))
        if not self.clinical and not self.fastq_dir:
            """Assert that the first FASTQ entry is non-empty and has the right format."""
            self._validate_fastq_format(row["fastq_1"])

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        print(filename.suffixes)
        assert any("".join(filename.suffixes) == extension for extension in self.VALID_FORMATS), (
            "The FASTQ file: {file} has an unrecognized extension: {suffix}\n It should be one of: {valid_format}".format(file=filename , valid_format=", ".join(self.VALID_FORMATS), suffix="".join(filename.suffixes))
        )

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if len(row) >= 2 and not self.clinical:
            row[self._single_col] = False
            first_col_suffix = Path(row[1]).suffixes[-2:]
            second_col_suffix = Path(row[2]).suffixes[-2:]
            if first_col_suffix != second_col_suffix:
                raise AssertionError("FASTQ pairs must have the same file extensions.")
        else:
            row[self._single_col] = True

        """Change barcode column name to sample."""
        if self.clinical:
            row["sample"] = row.pop("barcode")
            #change row["sample"] to lowercase
            row["sample"] = row["sample"].lower()

    def _check_barcode_format(self, row):
        pattern = r"^barcode([0-9][0-9]|100)$"
        if not re.match(pattern, row["sample"]):
                raise ValueError("Invalid barcode format for sample {row}. Expected format: barcode01-barcode09 or barcode10-barcode100.".format(row=row['sample']))
        return True

    def _create_fastq_path(self, row):
        print(row["sample"])
        print(row["status"])
        if not self.fastq_dir:
            return None
        
        fastq_dir = Path(self.fastq_dir)
        if not fastq_dir.is_dir():
            raise FileNotFoundError("The FASTQ directory does not exist: {fastq_dir}".format(fastq_dir=fastq_dir))
        
        if self.clinical:
            print(row["sample"])
            print(row["status"])
            filename = row["sample"]
            if row["status"] == "discontinued":
                return None
        else:
            filename = row["filename"]
        
        pattern = "{filename}.*".format(filename=filename)
        matching_files = list(fastq_dir.glob(pattern))
        
        if len(matching_files) == 0:
            warnings.warn("No matching FASTQ file found in the directory: {fastq_dir} for file: {filename}".format(fastq_dir=fastq_dir, filename=filename))
            row["status"] = "discontinued"
            return None
        elif len(matching_files) > 1:
            raise FileNotFoundError("Multiple matching FASTQ files found in the directory: {fastq_dir} for file: {filename}".format(fastq_dir=fastq_dir, filename=filename))
        
        fastq_file = matching_files[0]

        row["fastq_1"] = fastq_file
        self._validate_fastq_format(row["fastq_1"])
        
        return row

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename all samples to have a suffix of _T{n}, where n is the
        number of times the same sample exist, but with different FASTQ files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and FASTQ must be unique.")
        seen = Counter()
        for row in self.modified:
            if self.clinical:
                sample = row["specimen number"]
                seen[sample] += 1
                row["specimen number"] = "{sample}_T{seen}".format(sample=sample, seen=seen[sample])
            else:
                sample = row["sample"]
                seen[sample] += 1
                row["sample"] = "{sample}_T{seen}".format(sample=sample, seen=seen[sample])


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out, fastq_dir=None, clinical=False):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. The required columns depend on the conditional structure:
    - If the `clinical` flag is True, the samplesheet should contain columns "Specimen Number", "Barcode", and "Status".
    - If the `fastq_dir` is provided, the samplesheet should contain columns "Sample" and "FileName".
    - Otherwise, the samplesheet should contain columns "Sample" and "fastq_1".

    An additional column is added to record whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure, depending on the conditions:
        - If `clinical` is True:

            ```
            Specimen Number, Barcode, Status
            SAMPLE_1, BARCODE01, STATUS_1
            SAMPLE_2, BARCODE02, STATUS_2
            ```

        - If `fastq_dir` is provided:

            ```
            Sample, FileName
            SAMPLE_1, FILE_NAME_1
            SAMPLE_2, FILE_NAME_2
            ```

        - Otherwise:

            ```
            Sample, fastq_1
            SAMPLE_1, /path/to/SAMPLE_1.fastq[.gz]
            SAMPLE_2, /path/to/SAMPLE_2.fastq[.gz]
            ```

    """
    if clinical:
        required_columns = {"specimen number", "barcode", "status"}
    elif fastq_dir:
        required_columns = {"sample", "filename"}
    else:
        required_columns = {"sample", "fastq_1"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.

    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical("The sample sheet **must** contain these column headers: {req_cols}.".format(req_cols=req_cols))
            sys.exit(1)

        # Validate each row.
        checker = RowChecker(required_columns, fastq_dir=fastq_dir, clinical=clinical)
        for i, row in enumerate(reader):
            # if i == 0:
            #     print("Skipping the column names row.")
            #     continue
            try:
                print(row)
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical("{error} On line {i}.".format(error=str(error), i=i+2))
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    if clinical:
        header[header.index("barcode")] = "sample"
        header.append("fastq_1")
    
    header = ["sample", "fastq_1"] + [col for col in header if col not in ["sample", "fastq_1"]]
    
    header.insert(1, "single_end")
    header = [element.replace(" ", "_") for element in header]
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            new_row = {key.replace(' ', '_'): value for key, value in row.items()}
            writer.writerow(new_row)


def parse_args():
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "--fastq_dir",
        metavar="FASTQ_DIR",
        type=Path,
        help="Path to a directory with basecalled reads.",
    )
    parser.add_argument(
        "--clinical",
        action="store_true",
        help="Whether the samplesheet is for a clinical study.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args()


def convert_excel_to_csv(excel_file, clinical=False):
    try:
        # Validate file format
        excel_formats = pd.ExcelFile(excel_file).sheet_names
        if not excel_formats:
            logger.error("The input file does not contain any Excel sheets but the file has an xlsx extension. Check the file isn't corrupted.")
            sys.exit(1)

        # Perform conversion
        df = pd.read_excel(excel_file)
        df.columns = df.columns.str.lower()

        if clinical:
            df = df.iloc[:, :14]
            df_cleaned = df.dropna(how='all')

        csv_file = excel_file.with_suffix(".csv")
        df_cleaned.to_csv(csv_file, index=False)
        return csv_file
    except Exception as e:
        logger.error("An error occurred during the Excel to CSV conversion process:")
        logger.error(str(e))
        sys.exit(1)


def main():
    """Coordinate argument parsing and program execution."""
    args = parse_args()
    print(args)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error("The given input file {file} was not found!".format(file=args.file_in))
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)

    # Convert Excel file to CSV format
    if args.file_in.suffix == ".xlsx":
        csv_file = convert_excel_to_csv(args.file_in, args.clinical)
        args.file_in = csv_file

    print(args.clinical)   
    check_samplesheet(args.file_in, args.file_out, args.fastq_dir, args.clinical)


if __name__ == "__main__":
    sys.exit(main())
