#!/usr/bin/env python3

import sys
import argparse
from pathlib import Path

# Add scripts directory to path
PROJECT_DIR = Path.home() / "resurrection_lncrna_pipeline"
sys.path.insert(0, str(PROJECT_DIR / "scripts"))

from scripts.data_loader import DataLoader
from scripts.analysis import main as run_analysis
from scripts.visualization import Visualizer

def main():
    """Main entry point for the pipeline"""
    
    parser = argparse.ArgumentParser(
        description="Resurrection Plant lncRNA Discovery Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run.py --check-structure
  python run.py --load-data
  python run.py --analyze
  python run.py --visualize
  python run.py --full
        """
    )
    
    parser.add_argument('--check-structure', action='store_true',
                       help='Check project directory structure')
    parser.add_argument('--load-data', action='store_true',
                       help='Load and display available data files')
    parser.add_argument('--analyze', action='store_true',
                       help='Run differential expression analysis')
    parser.add_argument('--visualize', action='store_true',
                       help='Create visualizations (requires analysis results)')
    parser.add_argument('--full', action='store_true',
                       help='Run full pipeline (check + analyze + visualize)')
    
    args = parser.parse_args()
    
    # If no arguments, show help and status
    if not any(vars(args).values()):
        parser.print_help()
        print("\n" + "="*80)
        print("QUICK STATUS")
        print("="*80)
        loader = DataLoader()
        loader.check_project_structure()
        loader.list_results()
        return 0
    
    # Check structure
    if args.check_structure or args.full:
        print("="*80)
        print("CHECKING PROJECT STRUCTURE")
        print("="*80)
        print()
        loader = DataLoader()
        loader.check_project_structure()
        loader.list_results()
        print()
    
    # Load data
    if args.load_data or args.full:
        print("="*80)
        print("LOADING DATA")
        print("="*80)
        print()
        loader = DataLoader()
        
        tpm = loader.load_tpm_matrix()
        coding = loader.load_coding_potential()
        annotation = loader.load_annotation()
        
        print()
    
    # Analyze
    if args.analyze or args.full:
        print("="*80)
        print("RUNNING ANALYSIS")
        print("="*80)
        print()
        exit_code = run_analysis()
        if exit_code != 0:
            print("Analysis failed!")
            return exit_code
        print()
    
    # Visualize
    if args.visualize or args.full:
        print("="*80)
        print("CREATING VISUALIZATIONS")
        print("="*80)
        print()
        
        try:
            loader = DataLoader()
            tpm = loader.load_tpm_matrix()
            de_results = loader.load_differential_expression()
            sig_lncrnas = loader.load_significant_lncrnas()
            
            if de_results is not None and sig_lncrnas is not None:
                visualizer = Visualizer()
                
                if tpm is not None:
                    visualizer.volcano_plot(de_results)
                    visualizer.expression_boxplot(tpm, sig_lncrnas)
                
                visualizer.summary_statistics(de_results)
                
                print("\nVisualizations created successfully!")
            else:
                print("Could not load results. Please run --analyze first.")
                return 1
        
        except Exception as e:
            print(f"Visualization failed: {e}")
            return 1
        
        print()
    
    print("="*80)
    print("PIPELINE COMPLETE")
    print("="*80)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
