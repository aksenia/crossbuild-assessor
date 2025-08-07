#!/usr/bin/env python3
"""
VCF Liftover: Tool to annotate VCF with coordinates from another reference genome.

This script annotates a VCF file with coordinates from another reference genome 
(e.g., hg38) using CrossMap for coordinate conversion.
"""

import argparse
import subprocess
import tempfile
from pathlib import Path
import pysam
import logging
import shutil
import re
import time
from typing import NamedTuple, Generator, List, Dict, Tuple, Optional, Union
import os
from dataclasses import dataclass

# Configuration object for genome-specific parameters
@dataclass
class GenomeConfig:
    """Configuration for genome-specific parameters."""
    use_chr_prefix: bool = True
    default_chr_prefix: str = "chr"
    
    def format_chrom(self, chrom: str) -> str:
        """Format chromosome name according to configuration."""
        if self.use_chr_prefix and not chrom.startswith(self.default_chr_prefix):
            return f"{self.default_chr_prefix}{chrom}"
        elif not self.use_chr_prefix and chrom.startswith(self.default_chr_prefix):
            return chrom[len(self.default_chr_prefix):]
        return chrom
    
    def normalize_chrom(self, chrom: str) -> str:
        """Normalize chromosome name for consistent comparison."""
        if chrom.startswith(self.default_chr_prefix):
            return chrom[len(self.default_chr_prefix):]
        return chrom


class BedLine(NamedTuple):
    """Represents a line in a BED file."""
    chrom: str
    start: int
    end: int
    name: str

    def __str__(self) -> str:
        """Convert to string representation for BED file."""
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}"


class LiftedBedLine(NamedTuple):
    """Represents a lifted coordinate in BED format."""
    orig_chrom: str
    orig_start: int
    orig_end: int
    hg38_chrom: str
    hg38_start: Union[int, str]  # Allow string for "."
    hg38_end: Union[int, str]    # Allow string for "."
    status: str
    
    def __str__(self) -> str:
        """Convert to string representation for BED file."""
        # For FAILED records, use "." for all hg38 fields
        if self.status == "FAILED":
            return (f"{self.orig_chrom}\t{self.orig_start}\t{self.orig_end}\t"
                    f".\t.\t.\t.\t{self.status}")
        else:
            # Normal formatting for successful mappings
            hg38_coord = f"{self.hg38_chrom}:{self.hg38_start+1 if isinstance(self.hg38_start, int) else '.'}-{self.hg38_end}"
            return (f"{self.orig_chrom}\t{self.orig_start}\t{self.orig_end}\t"
                    f"{hg38_coord}\t"
                    f"{self.hg38_chrom}\t{self.hg38_start}\t{self.hg38_end}\t{self.status}")

class VcfAnnotation(NamedTuple):
    """Represents an annotated VCF entry."""
    fields: List[str]
    hg38_annotations: Dict[str, str]

    def to_line(self) -> str:
        """Convert to string representation for VCF file."""
        fields = self.fields.copy()
        info = fields[7]

        # Append annotations to INFO
        if self.hg38_annotations:
            new_info = [f"{key}={value}" for key, value in self.hg38_annotations.items()]
            updated_info = f"{info};{';'.join(new_info)}" if info else ';'.join(new_info)
        else:
            updated_info = info

        fields[7] = updated_info
        return '\t'.join(fields)

class CommandRunner:
    """Handles running external commands with error handling and retries."""
    
    def __init__(self, max_retries: int = 3, retry_delay: int = 2):
        self.max_retries = max_retries
        self.retry_delay = retry_delay
    
    def run(self, cmd: List[str], check: bool = True, capture_output: bool = True) -> subprocess.CompletedProcess:
        """
        Run a command with retries and exponential backoff.
        
        Args:
            cmd: Command to run as a list of strings
            check: Whether to raise an exception on failure
            capture_output: Whether to capture stdout and stderr
            
        Returns:
            CompletedProcess instance with output
            
        Raises:
            subprocess.CalledProcessError: If command fails and check is True
        """
        for attempt in range(self.max_retries):
            try:
                result = subprocess.run(
                    cmd, 
                    check=check, 
                    capture_output=capture_output, 
                    text=True
                )
                return result
            except subprocess.CalledProcessError as e:
                if attempt < self.max_retries - 1:
                    delay = self.retry_delay * (2 ** attempt)  # Exponential backoff
                    logging.warning(
                        f"Command failed (attempt {attempt+1}/{self.max_retries}): {' '.join(cmd)}"
                        f"\nRetrying in {delay} seconds..."
                    )
                    time.sleep(delay)
                else:
                    logging.error(
                        f"Command failed after {self.max_retries} attempts: {' '.join(cmd)}"
                        f"\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}"
                    )
                    raise


class FileHandler:
    """Handles file operations with safety features."""
    
    @staticmethod
    def atomic_write(path: Path, content: str, mode: str = 'w') -> None:
        """Write content to a file atomically."""
        temp_path = path.with_suffix('.tmp')
        try:
            with open(temp_path, mode) as f:
                f.write(content)
            # Atomic rename
            temp_path.replace(path)
        except Exception as e:
            logging.error(f"Error writing file {path}: {e}")
            if temp_path.exists():
                temp_path.unlink()
            raise
    
    @staticmethod
    def atomic_copy(src: Path, dst: Path) -> None:
        """Copy a file atomically."""
        temp_dst = dst.with_suffix('.tmp')
        try:
            shutil.copy(str(src), str(temp_dst))
            # Atomic rename
            temp_dst.replace(dst)
        except Exception as e:
            logging.error(f"Error copying {src} to {dst}: {e}")
            if temp_dst.exists():
                temp_dst.unlink()
            raise
    
    @staticmethod
    def ensure_bgzipped(file_path: Path) -> Path:
        """
        Ensure a file is bgzip-compressed. If not, compress it.
        
        Args:
            file_path: Path to the file
            
        Returns:
            Path to the compressed file
        """
        if FileHandler.is_bgzf_compressed(file_path):
            return file_path
        
        compressed_path = file_path.with_suffix(file_path.suffix + '.gz')
        cmd_runner = CommandRunner()
        cmd_runner.run(['bgzip', '-c', str(file_path)], check=True)
        
        if not compressed_path.exists():
            raise FileNotFoundError(f"Failed to compress {file_path}")
        
        return compressed_path
    
    @staticmethod
    def is_bgzf_compressed(file_path: Path) -> bool:
        """
        Check if a file is BGZF compressed.
        
        Args:
            file_path: Path to the file
            
        Returns:
            True if the file is BGZF compressed, False otherwise
        """
        try:
            # First, check if the file has the gzip magic bytes
            with open(str(file_path), 'rb') as f:
                magic = f.read(2)
                if magic != b'\x1f\x8b':  # Gzip magic bytes
                    return False
                    
            # If it has the right magic bytes, try pysam's test
            with pysam.BGZFile(str(file_path)) as _:
                return True
        except (OSError, IOError):
            return False

    @staticmethod
    def validate_file_exists(file_path: Path) -> None:
        """
        Validate that a file exists.
        
        Args:
            file_path: Path to the file
            
        Raises:
            FileNotFoundError: If the file does not exist
        """
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
    
    @staticmethod
    def validate_is_file(file_path: Path) -> None:
        """
        Validate that a path is a file.
        
        Args:
            file_path: Path to validate
            
        Raises:
            ValueError: If the path is not a file
        """
        if not file_path.is_file():
            raise ValueError(f"Not a file: {file_path}")

    @staticmethod
    def debug_bed_file(bed_path: Path, max_lines: int = 5, label: str = "BED file", log_level: int = logging.DEBUG) -> None:
        """
        Print contents of a BED file for debugging in a clean format.
        
        Args:
            bed_path: Path to the BED file
            max_lines: Maximum number of lines to print
            label: Label to use in the log messages
            log_level: Logging level to use (default: DEBUG)
        """
        try:
            if not os.path.exists(bed_path):
                logging.debug(f"{label} does not exist: {bed_path}")
                return
                
            if os.path.getsize(bed_path) == 0:
                logging.debug(f"{label} is empty: {bed_path}")
                return
                
            line_count = 0
            with open(bed_path, 'r') as f:
                lines = []
                for i, line in enumerate(f):
                    if i >= max_lines:
                        break
                    lines.append(line.strip())
                    line_count += 1
                
                total_lines = line_count
                if i == max_lines - 1:
                    # We hit the limit, count the rest
                    for _ in f:
                        total_lines += 1
                
                logging.debug(f"{label} {bed_path} contents ({len(lines)} of {total_lines} lines):")
                for i, line in enumerate(lines):
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        logging.debug(f"  [{i+1}] {parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}")
                    else:
                        logging.debug(f"  [{i+1}] {line}")
        except Exception as e:
            logging.error(f"Error reading {label}: {e}")

class BEDGenerator:
    """Handles conversion from VCF to BED format."""
    
    def __init__(self, genome_config: GenomeConfig):
        self.genome_config = genome_config
    

    def vcf_to_bed(self, input_vcf: Path) -> Generator[BedLine, None, None]:
        """
        Generate BED lines from VCF records with proper handling of structural variants.

        Args:
            input_vcf: Path to the input VCF file

        Returns:
            Generator of BedLine objects
        """
        logging.info(f"Processing VCF to generate BED data: {input_vcf}")

        # Track unique entries to avoid duplicates
        seen_entries = set()

        try:
            with pysam.VariantFile(str(input_vcf)) as vcf:
                line_count = 0
                duplicate_count = 0

                for record in vcf:
                    chrom = record.chrom
                    formatted_chrom = self.genome_config.format_chrom(chrom)
                    start = record.pos - 1  # Convert 1-based VCF to 0-based BED

                    # Get END from INFO - critical for structural variants
                    end = None
                    if 'END' in record.info:
                        end = record.info['END']
                    elif 'SVLEN' in record.info:
                        # If we have SVLEN but not END, calculate END
                        svlen = record.info['SVLEN']
                        if isinstance(svlen, tuple):
                            svlen = svlen[0]
                        # Control the Type of SVLEN - String or int    
                        if isinstance(svlen, str):
                            svlen = int(svlen)
    
                        # SVLEN is negative for deletions, positive for insertions
                        if svlen < 0:
                            end = record.pos - svlen
                        else:
                            end = record.pos + svlen

                    # If no END info, use position + ref length - 1
                    if end is None:
                        end = record.pos + len(record.ref) - 1

                    # Create a unique key for this entry
                    entry_key = f"{formatted_chrom}:{start}:{end}:{chrom}:{record.pos}:{end}"

                    # Skip if we've seen this entry before
                    if entry_key in seen_entries:
                        duplicate_count += 1
                        if logging.getLogger().isEnabledFor(logging.DEBUG) and duplicate_count <= 5:
                            logging.debug(f"Skipping duplicate entry: {formatted_chrom}\t{start}\t{end}\t{chrom}:{record.pos}:{end}")
                        continue

                    seen_entries.add(entry_key)
                    name = f"{chrom}:{record.pos}:{end}"

                    yield BedLine(formatted_chrom, start, end, name)

                    # Log for debug
                    line_count += 1
                    if logging.getLogger().isEnabledFor(logging.DEBUG) and line_count <= 10:
                        logging.debug(f"Generated BED line: {formatted_chrom}\t{start}\t{end}\t{name}")

            logging.info(f"Generated {line_count} unique BED lines")
        except Exception as e:
            logging.error(f"Error processing VCF file {input_vcf}: {e}")
            raise

class BEDSorter:
    """Handles sorting of BED files."""
    
    def __init__(self, cmd_runner: CommandRunner):
        self.cmd_runner = cmd_runner

    def sort_bed_file(self, bed_path: Path, output_path: Optional[Path] = None, compress: bool = False) -> Path:
        """Sort a BED file for tabix compatibility."""
        if not os.path.exists(bed_path):
            raise FileNotFoundError(f"BED file not found: {bed_path}")
            
        if output_path is None:
            output_path = bed_path
       
        # Handle empty input files
        if os.path.getsize(bed_path) == 0:
            logging.warning(f"Input BED file is empty: {bed_path}")
            if output_path is not None and output_path != bed_path:
                # Create an empty output file
                with open(output_path, 'w') as f:
                    pass
                return output_path
            return bed_path

        try:
            # Use subprocess with the correct sort parameters
            temp_sorted_path = output_path.with_suffix('.sorted.tmp')
            
            # Run the sort command directly
            self.cmd_runner.run([
                'sort', 
                '-k1,1V',  # Version sort for chromosome
                '-k2,2n',  # Numeric sort for position
                '-o', str(temp_sorted_path),
                str(bed_path)
            ])
            
            # Verify the sorted file exists and has content
            if not os.path.exists(temp_sorted_path) or os.path.getsize(temp_sorted_path) == 0:
                raise RuntimeError(f"Sort produced empty output file: {temp_sorted_path}")
            
            # Move the sorted file to the output path
            shutil.move(str(temp_sorted_path), str(output_path))
           
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                FileHandler.debug_bed_file(bed_path=output_path, max_lines=5, label="Sorted BED file", log_level=logging.DEBUG)

            # Compress if requested
            if compress:
                self.cmd_runner.run(['bgzip', '-f', str(output_path)])
                output_path = Path(str(output_path) + '.gz')
            
            logging.info(f"Successfully sorted BED file: {output_path}")
            return output_path
            
        except Exception as e:
            logging.error(f"Error sorting BED file {bed_path}: {e}")
            # Clean up any temporary files
            temp_sorted_path = output_path.with_suffix('.sorted.tmp')
            if os.path.exists(temp_sorted_path):
                os.unlink(temp_sorted_path)
            raise

class LiftOverHandler:
    """Handles coordinate lifting using CrossMap."""
    
    def __init__(self, chain_file: Path, cmd_runner: CommandRunner, tmp_dir: Path):
        self.chain_file = chain_file
        self.cmd_runner = cmd_runner
        self.tmp_dir = tmp_dir
        
        # Validate chain file
        FileHandler.validate_file_exists(chain_file)
        FileHandler.validate_is_file(chain_file)
    
    def lift_coordinates(self, bed_path: Path) -> Tuple[Path, Path]:
        """
        Run CrossMap to lift coordinates.
        
        Args:
            bed_path: Path to the BED file
            
        Returns:
            Tuple of (lifted_path, unlifted_path)
        """
        logging.info(f"Running CrossMap on BED file: {bed_path}")
        
        lifted_path = self.tmp_dir / 'lifted.bed'
        unlifted_path = self.tmp_dir / 'unlifted.bed'

        # Count lines to verify content
        line_count = 0
        with open(bed_path, 'r') as f:
            for line in f:
                if line.strip():  # Skip empty lines
                    line_count += 1
        
        logging.info(f"Input BED file contains {line_count} lines")
        
        if line_count == 0:
            logging.error(f"BED file is empty: {bed_path}")
            with open(lifted_path, 'w') as f:
                pass
            with open(unlifted_path, 'w') as f:
                pass
            return lifted_path, unlifted_path

        # Run CrossMap
        try:
            result = self.cmd_runner.run([
                'CrossMap', 'bed',
                '--unmap-file', str(unlifted_path),
                str(self.chain_file),
                str(bed_path),
                str(lifted_path)
            ], capture_output=True)
            
            # Only log detailed information in debug mode
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                # Log CrossMap errors
                if result.stderr:
                    logging.debug(f"CrossMap stderr: {result.stderr}")

                # Check lifted file
                FileHandler.debug_bed_file(lifted_path, max_lines=5, label="Lifted coordinates", log_level=logging.DEBUG)

                # Check unlifted file
                FileHandler.debug_bed_file(unlifted_path, max_lines=5, label="Unlifted coordinates", log_level=logging.DEBUG)
                      
            # Verify CrossMap output in all modes
            if not os.path.exists(lifted_path):
                logging.error("CrossMap output file does not exist")
        except Exception as e:
            logging.error(f"CrossMap failed: {e}")
            with open(lifted_path, 'w') as f:
                pass
        
        return lifted_path, unlifted_path


    def run_region_mapping(self, multi_path: Path) -> Path:
            """
            Run CrossMap region mapping for multi-mapped coordinates.
            
            Args:
                multi_path: Path to the BED file with multi-mapped coordinates
                
            Returns:
                Path to the region-mapped BED file
            """
            region_path = self.tmp_dir / 'region_mapped.bed'
            
            self.cmd_runner.run([
                'CrossMap', 'region',
                str(self.chain_file),
                str(multi_path),
                str(region_path)
            ])
            
            return region_path


class CoordinateProcessor:
    """Processes lifted coordinates for annotation."""
    
    def __init__(self, tmp_dir: Path, genome_config: GenomeConfig):
        self.tmp_dir = tmp_dir
        self.genome_config = genome_config
    
    def process_lifted_coordinates(self, 
                                  lifted_path: Path,
                                  unlifted_path: Path,
                                  original_bed_path: Path,
                                  cmd_runner: CommandRunner,
                                  chain_file: Path) -> Generator[LiftedBedLine, None, None]:
        """
        Process lifted coordinates and generate annotated bed lines.
        
        Args:
            lifted_path: Path to the lifted BED file
            original_bed_path: Path to the original BED file
            cmd_runner: CommandRunner instance
            chain_file: Path to the chain file
            
        Returns:
            Generator of LiftedBedLine objects
        """
        # First pass: Load and analyze lifted coordinates
        logging.info("Analyzing lifted coordinates for multi-mappings")
        
        # Load all lifted coordinates - use streaming for larger files
        mapping_counts = {}
        lifted_entries = []
        
        with open(lifted_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    entry = BedLine(parts[0], int(parts[1]), int(parts[2]), parts[3])
                    lifted_entries.append(entry)
                    
                    # Count mappings
                    original_coords = entry.name
                    mapping_counts[original_coords] = mapping_counts.get(original_coords, 0) + 1
        
        # Identify multi-mapped coordinates
        multi_coords = {coord for coord, count in mapping_counts.items() if count > 1}
        
        # Stats
        total_entries = len(lifted_entries)
        unique_sources = len(mapping_counts)
        multi_mapped = len(multi_coords)
        
        logging.info(f"Total entries in lifted BED: {total_entries}")
        logging.info(f"Unique source coordinates: {unique_sources}")
        logging.info(f"Multi-mapped source coordinates: {multi_mapped}")
       
        # Read unlifted coordinates
        unlifted_coords = set()
        if os.path.exists(unlifted_path) and os.path.getsize(unlifted_path) > 0:
            with open(unlifted_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        unlifted_coords.add(parts[3])  # Add name to unlifted set
            
            logging.info(f"Found {len(unlifted_coords)} unlifted coordinates")

        # Process unique mappings
        processed_coords = set()
        
        for entry in lifted_entries:
            original_coords = entry.name
            
            if mapping_counts[original_coords] == 1:
                # Handle unique mapping
                hg38_chr, hg38_start, hg38_end = entry.chrom, entry.start, entry.end
                original_chr, original_start, original_end = original_coords.split(':')
    
                # If we need to remove chr prefix
                if not self.genome_config.use_chr_prefix:
                    hg38_chr = self.genome_config.normalize_chrom(hg38_chr)

                yield LiftedBedLine(
                    orig_chrom=original_chr,
                    orig_start=int(original_start)-1,
                    orig_end=int(original_end),
                    hg38_chrom=hg38_chr,
                    hg38_start=hg38_start,
                    hg38_end=hg38_end,
                    status="UNIQUE"
                )
                
                processed_coords.add(original_coords)
        
        # Handle multi-mapped coordinates using region mapping
        if multi_coords:
            # Create temporary BED file for multi-mapped coordinates
            multi_path = self.tmp_dir / 'multi_mappings_original.bed'
            with open(multi_path, 'w') as f:
                with open(original_bed_path, 'r') as original_file:
                    for line in original_file:
                        parts = line.strip().split('\t')
                        if len(parts) == 4 and parts[3] in multi_coords:
                            f.write(line)
            
            # Run CrossMap region
            liftover_handler = LiftOverHandler(chain_file, cmd_runner, self.tmp_dir)
            region_path = liftover_handler.run_region_mapping(multi_path)
            
            # Process region mappings
            region_processed = set()
            with open(region_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    
                    if len(parts) < 5:
                        continue
                    
                    hg38_chr, hg38_start, hg38_end, original_coords = parts[0], parts[1], parts[2], parts[3]
                    
                    if original_coords in processed_coords:
                        continue
                    
                    processed_coords.add(original_coords)
                    region_processed.add(original_coords)
                    original_chr, original_start, original_end = original_coords.split(':')
                    
                    yield LiftedBedLine(
                        orig_chrom=original_chr,
                        orig_start=int(original_start)-1,
                        orig_end=int(original_end),
                        hg38_chrom=hg38_chr,
                        hg38_start=int(hg38_start),
                        hg38_end=int(hg38_end),
                        status="REGION"
                    )
            
            # Handle failed mappings
            missing_coords = multi_coords - region_processed
            for original_coords in missing_coords:
                if original_coords in processed_coords:
                    continue
                
                original_chr, original_start, original_end = original_coords.split(':')
                
                yield LiftedBedLine(
                    orig_chrom=original_chr,
                    orig_start=int(original_start)-1,
                    orig_end=int(original_end),
                    hg38_chrom=".",
                    hg38_start=0,
                    hg38_end=0,
                    status="FAILED"
                )
                
                processed_coords.add(original_coords)

        # Process unlifted coordinates
        for coord in unlifted_coords:
            if coord in processed_coords:
                continue

            try:
                original_chr, original_start, original_end = coord.split(':')

                yield LiftedBedLine(
                    orig_chrom=original_chr,
                    orig_start=int(original_start)-1,
                    orig_end=int(original_end),
                    hg38_chrom=".",
                    hg38_start=".",
                    hg38_end=".",
                    status="FAILED"
                )

                processed_coords.add(coord)
            except Exception as e:
                logging.error(f"Error processing unlifted coordinate {coord}: {e}")
            
        logging.info(f"Processed {len(processed_coords)} coordinates out of {total_entries} lifted entries")
        if not processed_coords:
            logging.warning("No coordinates were processed! Check the lifted file content.")
        # Debug: inspect the lifted file
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            if os.path.exists(lifted_path) and os.path.getsize(lifted_path) > 0:
                FileHandler.debug_bed_file(lifted_path, max_lines=5, label="Sorted BED file", log_level=logging.DEBUG)

class VCFHeaderProcessor:
    """Handles VCF header processing."""
    
    @staticmethod
    def get_vcf_header(input_vcf: Path) -> List[str]:
        """
        Extract header lines from a VCF file.
        
        Args:
            input_vcf: Path to the VCF file
            
        Returns:
            List of header lines
        """
        header_lines = []
        
        # Read original header
        is_compressed = FileHandler.is_bgzf_compressed(input_vcf)
        
        if is_compressed:
            with pysam.BGZFile(str(input_vcf), 'r') as f:
                for line in f:
                    line = line.decode('utf-8').strip()
                    header_lines.append(line)
                    if line.startswith('#CHROM'):
                        break
        else:
            with open(input_vcf, 'r') as f:
                for line in f:
                    line = line.strip()
                    header_lines.append(line)
                    if line.startswith('#CHROM'):
                        break
        
        return header_lines
    
    @staticmethod
    def modify_header(header_lines: List[str], new_info_fields: List[str]) -> str:
        """
        Modify a VCF header with new INFO fields.
        
        Args:
            header_lines: List of header lines
            new_info_fields: List of new INFO field definitions
            
        Returns:
            Modified header as a string
        """
        # Find #CHROM line index
        try:
            chrom_idx = next(i for i, line in enumerate(header_lines) if line.startswith('#CHROM'))
        except StopIteration:
            raise ValueError("VCF header is missing #CHROM line")
        
        # Insert new INFO fields before #CHROM
        modified_header = header_lines.copy()
        for info in reversed(new_info_fields):
            modified_header.insert(chrom_idx, info)
        
        # Join with newlines
        return '\n'.join(modified_header)


class VCFAnnotator:
    """Handles annotation of VCF files with lifted coordinates."""
    
    def __init__(self, genome_config: GenomeConfig):
        self.genome_config = genome_config
    
    def get_vcf_annotations(self, input_vcf: Path, annotated_bed: Path) -> Generator[VcfAnnotation, None, None]:
        """
        Generate annotated VCF lines.
        
        Args:
            input_vcf: Path to the input VCF file
            annotated_bed: Path to the annotated BED file
            
        Returns:
            Generator of VcfAnnotation objects
        """
        logging.info(f"Starting VCF annotation with BED file: {annotated_bed}")
        annotation_count = 0

        is_compressed = FileHandler.is_bgzf_compressed(input_vcf)
        
        try:
            with pysam.TabixFile(str(annotated_bed)) as bed:
                if is_compressed:
                    with pysam.BGZFile(str(input_vcf), 'r') as f:
                        # Skip header lines
                        for line in f:
                            line = line.decode('utf-8').strip()
                            if line.startswith('#CHROM'):
                                break
                        
                        # Process data lines
                        for line in f:
                            line = line.decode('utf-8').strip()
                            yield self._generate_vcf_annotation(line, bed)
                else:
                    with open(input_vcf, 'r') as f:
                        # Skip header lines
                        for line in f:
                            line = line.strip()
                            if line.startswith('#CHROM'):
                                break
                        
                        # Process data lines
                        for line in f:
                            line = line.strip()
                            yield self._generate_vcf_annotation(line, bed)
        except Exception as e:
            logging.error(f"Error getting VCF annotations: {e}")
            raise
    
    def _generate_vcf_annotation(self, line: str, bed: pysam.TabixFile) -> VcfAnnotation:
        """
        Generate a VCF annotation from a line.
        
        Args:
            line: VCF line
            bed: TabixFile object for the annotated BED file
            
        Returns:
            VcfAnnotation object
        """
        fields = line.split('\t')
        if len(fields) < 8:
            return VcfAnnotation(fields, {})
        
        # Extract fields
        chrom, pos, id_val, ref, alt, qual, filter_val, info = fields[:8]
        pos_int = int(pos)
        
        # Get END value from INFO
        end = pos_int
        for item in info.split(';'):
            if item.startswith('END='):
                try:
                    end = int(item.split('=')[1])
                except (ValueError, IndexError):
                    pass
                break
        
        # Get hg38 annotations
        hg38_annotations = {}
        try:
            # Convert to 0-based for BED
            region_start = pos_int - 1
            
            # Try both with and without chr prefix
            bed_entries = []
            try:
                bed_entries = list(bed.fetch(chrom, region_start, end))
            except ValueError:
                # Try with chr prefix
                try:
                    bed_entries = list(bed.fetch(self.genome_config.format_chrom(chrom), region_start, end))
                except ValueError:
                    # Try without chr prefix
                    try:
                        plain_chrom = self.genome_config.normalize_chrom(chrom)
                        bed_entries = list(bed.fetch(plain_chrom, region_start, end))
                    except ValueError:
                        pass
            
            if bed_entries:
                bed_entry = bed_entries[0].split('\t')
        
                # Mapping status is in the last field
                mapping_status = bed_entry[7] if len(bed_entry) > 7 else "UNIQUE"
                
                # Handle FAILED entries differently
                if mapping_status == "FAILED":
                    hg38_annotations = {
                        'hg38_map': "FAILED"
                    }
                else:
                    # Strip "chr" prefix using regex
                    hg38_chr = bed_entry[4]
                    if not self.genome_config.use_chr_prefix:
                        hg38_chr = self.genome_config.normalize_chrom(hg38_chr)


                    # Handle hg38_coord field - strip "chr" prefix
                    
                    hg38_coord = bed_entry[3]
                    if not self.genome_config.use_chr_prefix:
                        # Extract and format the coordinate part
                        if ':' in hg38_coord:
                            coord_parts = hg38_coord.split(':')
                            chrom_part = self.genome_config.normalize_chrom(coord_parts[0])
                            hg38_coord = f"{chrom_part}:{coord_parts[1]}"

                    hg38_annotations = {
                        'hg38_chr': hg38_chr,
                        'hg38_start': str(int(bed_entry[5]) + 1) if bed_entry[5] != "." else ".",
                        'hg38_end': bed_entry[6],
                        'hg38_coord': hg38_coord,
                        'hg38_map': mapping_status
                    }

        except Exception as e:
            logging.warning(f"Error fetching hg38 annotations for {chrom}:{pos}: {e}")
        
        return VcfAnnotation(fields, hg38_annotations)


class Pipeline:
    """Main pipeline for VCF annotation."""
    
    def __init__(self, input_vcf: Path, chain_file: Path, outfile: Path, debug: bool, output_bed: bool = False):
        self.input_vcf = input_vcf
        self.chain_file = chain_file
        self.outfile = outfile
        self.debug = debug
        self.output_bed = output_bed
        self.genome_config = GenomeConfig(use_chr_prefix=False)
        
        # Create temporary directory
        self.tmpdir = tempfile.TemporaryDirectory(prefix='tempdir-')
        self.tmp_path = Path(self.tmpdir.name)
        
        # Initialize helpers
        self.cmd_runner = CommandRunner()
        self.bed_generator = BEDGenerator(self.genome_config)
        self.bed_sorter = BEDSorter(self.cmd_runner)
        self.liftover_handler = LiftOverHandler(self.chain_file, self.cmd_runner, self.tmp_path)
        self.coordinate_processor = CoordinateProcessor(self.tmp_path, self.genome_config)
        self.vcf_annotator = VCFAnnotator(self.genome_config)
        
        # Validate input files
        self._validate_inputs()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.tmpdir.cleanup()
    
    def _validate_inputs(self):
        """Validate input files."""
        FileHandler.validate_file_exists(self.input_vcf)
        FileHandler.validate_file_exists(self.chain_file)
        
        # Check if output directory exists
        outdir = self.outfile.parent
        if not outdir.exists():
            logging.info(f"Creating output directory: {outdir}")
            outdir.mkdir(parents=True, exist_ok=True)
   
    def run(self):
        """Run the pipeline"""
        # Step 1: Generate BED file from VCF
        bed_file = self.tmp_path / 'original.bed'
        bed_entries = 0
        
        with open(bed_file, "w") as f:
            for bed_line in self.bed_generator.vcf_to_bed(self.input_vcf):
                f.write(str(bed_line) + "\n")
                bed_entries += 1
        
        logging.info(f"Generated {bed_entries} BED entries from input VCF")
        
        # Check if we have any entries to process
        if bed_entries == 0:
            logging.warning("Input VCF produced no BED entries. Creating empty output files.")
            # Create empty output VCF
            with open(self.outfile, 'w') as f:
                # Copy the header from the input VCF
                with pysam.VariantFile(str(self.input_vcf)) as vcf:
                    f.write(str(vcf.header))
            
            # Create empty output BED if requested
            if self.output_bed:
                bed_outfile = Path(str(self.outfile) + '.bed.gz')
                with open(bed_outfile, 'w') as f:
                    pass  # Create empty file
                
            logging.info(f"Created empty output files: {self.outfile}")
            return self.outfile
        
        # Step 2: Sort BED file - ensure output is NOT compressed for CrossMap
        sorted_bed = self.tmp_path / 'sorted_uncompressed.bed'
        self.bed_sorter.sort_bed_file(bed_file, sorted_bed)

        # Step 3: Run liftover with uncompressed file
        lifted_path, unlifted_path = self.liftover_handler.lift_coordinates(sorted_bed)
        # Step 4: Process lifted coordinates
        annotated_bed = self.tmp_path / 'lifted_annotated.bed'
        with open(annotated_bed, "w") as f:
            for lifted_line in self.coordinate_processor.process_lifted_coordinates(
                lifted_path, unlifted_path, bed_file, self.cmd_runner, self.chain_file
            ):
                f.write(str(lifted_line) + "\n")
        
        # Step 5: Sort the annotated BED file for tabix indexing
        # This is critical - BED must be sorted by chromosome, then start position
        sorted_annotated_bed = self.tmp_path / 'sorted_annotated.bed'
        self.bed_sorter.sort_bed_file(annotated_bed, sorted_annotated_bed)

        # Step 6: Compress and index the sorted BED file
        annotated_bed_gz = sorted_annotated_bed.with_suffix('.bed.gz')
        self.cmd_runner.run(['bgzip', '-f', str(sorted_annotated_bed)])
        pysam.tabix_index(str(annotated_bed_gz), preset="bed", force=True)
        logging.info(f"Created tabix index for {annotated_bed_gz}")
        # Check if the file exists and has content
        if os.path.exists(annotated_bed_gz) and os.path.getsize(annotated_bed_gz) > 0:
            logging.info(f"Annotated BED file size: {os.path.getsize(annotated_bed_gz)} bytes")
        else:
            logging.error(f"Annotated BED file is empty or missing!")
                
        # Step 6: Prepare VCF header
        header_lines = VCFHeaderProcessor.get_vcf_header(self.input_vcf)
        new_info_fields = [
            '##INFO=<ID=hg38_chr,Number=1,Type=String,Description="Chromosome in hg38">',
            '##INFO=<ID=hg38_start,Number=1,Type=Integer,Description="Start position in hg38">',
            '##INFO=<ID=hg38_end,Number=1,Type=Integer,Description="End position in hg38">',
            '##INFO=<ID=hg38_coord,Number=.,Type=String,Description="Coordinates in hg38">',
            '##INFO=<ID=hg38_map,Number=1,Type=String,Description="Mapping status: UNIQUE, REGION, FAILED">'
        ]
        modified_header = VCFHeaderProcessor.modify_header(header_lines, new_info_fields)
        
        # Step 7: Annotate VCF
        logging.info(f"Annotating VCF using annotated BED file: {annotated_bed_gz}")
        with open(self.outfile, 'w') as outfile:

            # Write modified header
            outfile.write(modified_header)
            # Get all annotations
            annotations = list(self.vcf_annotator.get_vcf_annotations(self.input_vcf, annotated_bed_gz))
            if annotations:
                # Write header/data separator
                outfile.write('\n')
                # Write all records except the last one with newlines
                for annotation in annotations[:-1]:
                    outfile.write(annotation.to_line() + '\n')
                # Write the last record without a trailing newline
                outfile.write(annotations[-1].to_line())

        # Step 8: Save BED file if requested
        if self.output_bed:
            bed_outfile = Path(str(self.outfile) + '.bed.gz')
            logging.info(f"Saving annotated BED file to: {bed_outfile}")
            FileHandler.atomic_copy(annotated_bed_gz, bed_outfile)
            # Copy index too
            index_file = annotated_bed_gz.with_suffix('.bed.gz.tbi')
            if index_file.exists():
                index_outfile = bed_outfile.with_suffix('.bed.gz.tbi')
                FileHandler.atomic_copy(index_file, index_outfile)
                
        logging.info(f"Successfully created: {self.outfile}")
        return self.outfile

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description='Annotate VCF with hg38 coordinates using CrossMap',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('input_vcf', type=Path, help='Input VCF file')
    parser.add_argument('chain_file', type=Path, help='Path to liftOver chain file')
    parser.add_argument('outfile', type=Path, help='Output VCF file')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('--output-bed', action='store_true', help='Also output annotated BED file')

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.DEBUG if args.debug else logging.INFO
    )

    try:
        with Pipeline(
            input_vcf=args.input_vcf,
            chain_file=args.chain_file,
            outfile=args.outfile,
            debug=args.debug,
            output_bed=args.output_bed
        ) as pipeline:
            pipeline.run()
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        return 1
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
