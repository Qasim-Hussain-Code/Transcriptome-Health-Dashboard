"""
Interactive Visualization Module for RNA-Seq QC Dashboard

This module provides interactive Plotly-based visualizations for RNA-Seq
quality control metrics, replacing static Seaborn/Matplotlib plots with
hover-enabled, zoomable HTML dashboards.

Visualizations included:
- Library size distribution (histogram)
- Detected genes (boxplot/strip)
- PCA scatter plot
- Mitochondrial content bar chart
- Combined HTML dashboard
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import os
import logging
from typing import Optional, Dict, Any, List, TYPE_CHECKING

if TYPE_CHECKING:
    from .dataset import RNASeqDataset


# Vibrant color palette for scientific visualization
COLORS = {
    "primary": "#00D4FF",       # Bright Cyan
    "secondary": "#FF6B9D",     # Vibrant Pink/Magenta
    "tertiary": "#FFB347",      # Warm Orange
    "quaternary": "#9B59B6",    # Rich Purple
    "warning": "#FF4757",       # Bright Red
    "success": "#2ED573",       # Vibrant Green
    "background": "#0D1117",    # GitHub Dark
    "surface": "#161B22",       # Elevated surface
    "text": "#F0F6FC",          # Bright text
    "grid": "#30363D",          # Subtle grid
    "gradient_start": "#667EEA", # Blue-purple
    "gradient_end": "#764BA2"    # Purple
}

# Colorful gradient scales for different plot types
GRADIENT_LIBRARY = "Viridis"      # For library size
GRADIENT_GENES = "Plasma"         # For detected genes  
GRADIENT_PCA = "Turbo"            # For PCA
GRADIENT_MT = "RdYlGn_r"          # For MT content (red=bad, green=good)


def create_dark_theme():
    """Create a consistent dark theme for all plots."""
    return {
        "paper_bgcolor": COLORS["background"],
        "plot_bgcolor": COLORS["surface"],
        "font": {"color": COLORS["text"], "family": "Inter, system-ui, sans-serif"},
        "xaxis": {
            "gridcolor": COLORS["grid"],
            "linecolor": COLORS["grid"],
            "tickcolor": COLORS["text"]
        },
        "yaxis": {
            "gridcolor": COLORS["grid"],
            "linecolor": COLORS["grid"],
            "tickcolor": COLORS["text"]
        }
    }


def plot_library_sizes_interactive(
    qc_df: pd.DataFrame,
    threshold: int = 20_000_000,
    title: str = "Distribution of Sequencing Depth Across Samples"
) -> go.Figure:
    """
    Create an interactive histogram of library sizes.
    
    Features:
    - Hover shows Patient ID and exact library size
    - Red threshold line for minimum recommended reads
    - KDE density overlay
    
    Args:
        qc_df: DataFrame with 'Library_Size' column
        threshold: Minimum recommended library size (red line)
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    fig = go.Figure()
    
    # Main histogram with gradient color
    fig.add_trace(go.Histogram(
        x=qc_df["Library_Size"],
        nbinsx=30,
        name="Samples",
        marker=dict(
            color=qc_df["Library_Size"],
            colorscale="Viridis",
            line=dict(width=1, color=COLORS["text"])
        ),
        opacity=0.9,
        hovertemplate="<b>Library Size:</b> %{x:,.0f}<br><b>Count:</b> %{y}<extra></extra>"
    ))
    
    # Threshold line
    fig.add_vline(
        x=threshold,
        line_dash="dash",
        line_color=COLORS["warning"],
        line_width=2,
        annotation_text=f"Min Threshold ({threshold/1e6:.0f}M)",
        annotation_position="top right",
        annotation_font_color=COLORS["warning"]
    )
    
    # Calculate statistics for annotation
    mean_lib = qc_df["Library_Size"].mean()
    median_lib = qc_df["Library_Size"].median()
    below_threshold = (qc_df["Library_Size"] < threshold).sum()
    
    # Stats annotation
    stats_text = (
        f"<b>Summary Statistics</b><br>"
        f"Mean: {mean_lib/1e6:.1f}M<br>"
        f"Median: {median_lib/1e6:.1f}M<br>"
        f"Below Threshold: {below_threshold}"
    )
    
    fig.add_annotation(
        x=0.98, y=0.95,
        xref="paper", yref="paper",
        text=stats_text,
        showarrow=False,
        align="right",
        bgcolor=COLORS["surface"],
        bordercolor=COLORS["grid"],
        borderwidth=1,
        font=dict(size=11)
    )
    
    # Apply theme
    theme = create_dark_theme()
    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis_title="Library Size (Total Reads)",
        yaxis_title="Number of Samples",
        **theme,
        showlegend=False,
        hovermode="x unified"
    )
    
    # Format x-axis with millions
    fig.update_xaxes(tickformat=".2s")
    
    return fig


def plot_detected_genes_interactive(
    qc_df: pd.DataFrame,
    title: str = "Gene Detection Complexity per Sample"
) -> go.Figure:
    """
    Create an interactive boxplot with strip overlay for detected genes.
    
    Features:
    - Boxplot showing distribution summary
    - Individual points with hover showing Sample ID
    - Outliers highlighted automatically
    
    Args:
        qc_df: DataFrame with 'Detected_Genes' column and sample IDs as index
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    # Reset index to get sample IDs as column
    df = qc_df.reset_index()
    df.columns = ["Sample_ID"] + list(df.columns[1:])
    
    fig = go.Figure()
    
    # Box plot with vibrant color
    fig.add_trace(go.Box(
        y=df["Detected_Genes"],
        name="Distribution",
        marker_color=COLORS["secondary"],
        fillcolor="rgba(255, 107, 157, 0.3)",  # Translucent pink
        line=dict(color=COLORS["secondary"], width=2),
        boxmean="sd",
        boxpoints="outliers",
        jitter=0.3,
        hovertemplate="<b>Detected Genes:</b> %{y:,.0f}<extra></extra>"
    ))
    
    # Overlay strip chart with gradient colors based on gene count
    fig.add_trace(go.Scatter(
        y=df["Detected_Genes"],
        x=np.random.normal(0, 0.04, len(df)),  # Jitter
        mode="markers",
        name="Samples",
        marker=dict(
            color=df["Detected_Genes"],
            colorscale="Plasma",
            size=8,
            opacity=0.8,
            line=dict(width=1, color="white")
        ),
        text=df["Sample_ID"],
        hovertemplate="<b>Sample:</b> %{text}<br><b>Detected Genes:</b> %{y:,.0f}<extra></extra>"
    ))
    
    # Stats annotation
    q1 = df["Detected_Genes"].quantile(0.25)
    q3 = df["Detected_Genes"].quantile(0.75)
    iqr = q3 - q1
    median = df["Detected_Genes"].median()
    
    stats_text = (
        f"<b>Summary</b><br>"
        f"Median: {median:,.0f}<br>"
        f"IQR: {iqr:,.0f}<br>"
        f"Q1: {q1:,.0f} | Q3: {q3:,.0f}"
    )
    
    fig.add_annotation(
        x=0.98, y=0.95,
        xref="paper", yref="paper",
        text=stats_text,
        showarrow=False,
        align="right",
        bgcolor=COLORS["surface"],
        bordercolor=COLORS["grid"],
        borderwidth=1,
        font=dict(size=11)
    )
    
    # Apply theme
    theme = create_dark_theme()
    fig.update_layout(
        title=dict(text=title, x=0.5),
        yaxis_title="Number of Genes Detected (>0 counts)",
        **theme,
        showlegend=False
    )
    
    fig.update_xaxes(showticklabels=False, showgrid=False)
    
    return fig


def plot_pca_interactive(
    coordinates: np.ndarray,
    explained_variance: np.ndarray,
    sample_ids: List[str],
    metadata: Optional[pd.DataFrame] = None,
    color_by: Optional[str] = None,
    title: str = "Principal Component Analysis of Sample Expression Profiles"
) -> go.Figure:
    """
    Create an interactive PCA scatter plot.
    
    Features:
    - Hover shows Sample ID and PC coordinates
    - Optional coloring by metadata variable
    - Explained variance in axis labels
    
    Args:
        coordinates: PCA coordinates array (n_samples × n_components)
        explained_variance: Variance explained by each PC
        sample_ids: List of sample identifiers
        metadata: Optional DataFrame with sample metadata for coloring
        color_by: Column name in metadata to color points by
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    # Create DataFrame for plotting
    df = pd.DataFrame({
        "PC1": coordinates[:, 0],
        "PC2": coordinates[:, 1],
        "Sample_ID": sample_ids
    })
    
    # Add metadata if provided
    if metadata is not None and color_by and color_by in metadata.columns:
        df["Color"] = metadata[color_by].values
        fig = px.scatter(
            df, x="PC1", y="PC2",
            color="Color",
            hover_data=["Sample_ID"],
            title=title,
            color_discrete_sequence=px.colors.qualitative.Set2
        )
    else:
        # Color by PC1 value for gradient effect
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df["PC1"],
            y=df["PC2"],
            mode="markers",
            marker=dict(
                color=df["PC1"],
                colorscale="Turbo",
                size=12,
                opacity=0.85,
                line=dict(width=1, color="white"),
                showscale=True,
                colorbar=dict(title="PC1", thickness=15)
            ),
            text=df["Sample_ID"],
            hovertemplate=(
                "<b>Sample:</b> %{text}<br>"
                "<b>PC1:</b> %{x:.2f}<br>"
                "<b>PC2:</b> %{y:.2f}<extra></extra>"
            )
        ))
    
    # Apply theme
    theme = create_dark_theme()
    
    pc1_var = explained_variance[0] * 100
    pc2_var = explained_variance[1] * 100
    
    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis_title=f"PC1 ({pc1_var:.1f}% variance)",
        yaxis_title=f"PC2 ({pc2_var:.1f}% variance)",
        **theme,
        showlegend=True if color_by else False
    )
    
    # Add variance explained annotation
    total_var = (explained_variance[0] + explained_variance[1]) * 100
    fig.add_annotation(
        x=0.02, y=0.98,
        xref="paper", yref="paper",
        text=f"<b>Total Variance Explained: {total_var:.1f}%</b>",
        showarrow=False,
        align="left",
        bgcolor=COLORS["surface"],
        bordercolor=COLORS["grid"],
        borderwidth=1,
        font=dict(size=11)
    )
    
    return fig


def plot_mt_content_interactive(
    qc_df: pd.DataFrame,
    threshold: float = 20.0,
    title: str = "Mitochondrial Gene Expression as Percentage of Total Reads"
) -> go.Figure:
    """
    Create an interactive bar chart of mitochondrial percentage.
    
    Features:
    - Hover shows Sample ID and exact MT%
    - Samples exceeding threshold highlighted in red
    - Threshold line with annotation
    
    Args:
        qc_df: DataFrame with 'MT_Percentage' column
        threshold: Maximum acceptable MT percentage
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    # Reset index for sample IDs
    df = qc_df.reset_index()
    df.columns = ["Sample_ID"] + list(df.columns[1:])
    
    # Sort by MT percentage
    df = df.sort_values("MT_Percentage", ascending=False)
    
    # Use gradient based on MT value - green to red scale
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=list(range(len(df))),
        y=df["MT_Percentage"],
        marker=dict(
            color=df["MT_Percentage"],
            colorscale="RdYlGn_r",  # Red-Yellow-Green reversed (green=low, red=high)
            line=dict(width=1, color="white")
        ),
        text=df["Sample_ID"],
        hovertemplate=(
            "<b>Sample:</b> %{text}<br>"
            "<b>MT%:</b> %{y:.2f}%<extra></extra>"
        )
    ))
    
    # Threshold line
    fig.add_hline(
        y=threshold,
        line_dash="dash",
        line_color=COLORS["warning"],
        line_width=2,
        annotation_text=f"Threshold ({threshold}%)",
        annotation_position="top right",
        annotation_font_color=COLORS["warning"]
    )
    
    # Count failed samples
    n_failed = (df["MT_Percentage"] > threshold).sum()
    
    fig.add_annotation(
        x=0.98, y=0.95,
        xref="paper", yref="paper",
        text=f"<b>High MT% Samples:</b> {n_failed}",
        showarrow=False,
        align="right",
        bgcolor=COLORS["surface"],
        bordercolor=COLORS["grid"],
        borderwidth=1,
        font=dict(size=11)
    )
    
    # Apply theme
    theme = create_dark_theme()
    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis_title="Samples (sorted by MT%)",
        yaxis_title="Mitochondrial Content (%)",
        **theme,
        showlegend=False
    )
    
    fig.update_xaxes(showticklabels=False)
    
    return fig


def create_qc_summary_table(qc_df: pd.DataFrame) -> go.Figure:
    """
    Create a summary statistics table.
    
    Args:
        qc_df: DataFrame with QC metrics
        
    Returns:
        Plotly Figure with table
    """
    # Calculate summary stats
    stats = {
        "Metric": ["Library Size", "Detected Genes", "MT Percentage"],
        "Mean": [
            f"{qc_df['Library_Size'].mean():,.0f}",
            f"{qc_df['Detected_Genes'].mean():,.0f}",
            f"{qc_df['MT_Percentage'].mean():.2f}%"
        ],
        "Median": [
            f"{qc_df['Library_Size'].median():,.0f}",
            f"{qc_df['Detected_Genes'].median():,.0f}",
            f"{qc_df['MT_Percentage'].median():.2f}%"
        ],
        "Min": [
            f"{qc_df['Library_Size'].min():,.0f}",
            f"{qc_df['Detected_Genes'].min():,.0f}",
            f"{qc_df['MT_Percentage'].min():.2f}%"
        ],
        "Max": [
            f"{qc_df['Library_Size'].max():,.0f}",
            f"{qc_df['Detected_Genes'].max():,.0f}",
            f"{qc_df['MT_Percentage'].max():.2f}%"
        ]
    }
    
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=list(stats.keys()),
            fill_color=COLORS["surface"],
            font=dict(color=COLORS["text"], size=12),
            align="left",
            line_color=COLORS["grid"]
        ),
        cells=dict(
            values=list(stats.values()),
            fill_color=COLORS["background"],
            font=dict(color=COLORS["text"], size=11),
            align="left",
            line_color=COLORS["grid"]
        )
    )])
    
    fig.update_layout(
        title=dict(text="QC Summary Statistics", x=0.5),
        paper_bgcolor=COLORS["background"],
        font=dict(color=COLORS["text"])
    )
    
    return fig


def create_interactive_dashboard(
    dataset: 'RNASeqDataset',
    output_dir: str,
    lib_threshold: int = 20_000_000,
    mt_threshold: float = 20.0,
    filename: str = "dashboard.html"
) -> str:
    """
    Generate a comprehensive HTML dashboard with all QC plots.
    
    Creates a single HTML file containing:
    - Summary statistics table
    - Library size distribution
    - Detected genes boxplot
    - PCA scatter (if available)
    - Mitochondrial content bar chart
    
    Args:
        dataset: RNASeqDataset instance with QC metrics calculated
        output_dir: Directory to save the dashboard
        lib_threshold: Library size threshold for flagging
        mt_threshold: MT percentage threshold for flagging
        filename: Output filename
        
    Returns:
        Path to saved HTML file
    """
    logger = dataset.logger if hasattr(dataset, 'logger') else logging.getLogger(__name__)
    
    # Ensure QC metrics are calculated
    if dataset.qc_metrics is None:
        dataset.calculate_qc_metrics()
    
    qc_df = dataset.qc_metrics
    
    # Create individual plots with their own threshold lines
    fig_lib = plot_library_sizes_interactive(qc_df, threshold=lib_threshold)
    fig_genes = plot_detected_genes_interactive(qc_df)
    fig_mt = plot_mt_content_interactive(qc_df, threshold=mt_threshold)
    fig_summary = create_qc_summary_table(qc_df)
    
    # PCA plot if available
    has_pca = dataset.pca_results is not None
    if has_pca:
        fig_pca = plot_pca_interactive(
            dataset.pca_results["coordinates"],
            dataset.pca_results["explained_variance"],
            dataset.sample_ids
        )
    
    # Create combined HTML dashboard using div layout
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    
    # Convert figures to HTML divs
    summary_html = fig_summary.to_html(full_html=False, include_plotlyjs=False)
    lib_html = fig_lib.to_html(full_html=False, include_plotlyjs=False)
    genes_html = fig_genes.to_html(full_html=False, include_plotlyjs=False)
    mt_html = fig_mt.to_html(full_html=False, include_plotlyjs=False)
    pca_html = fig_pca.to_html(full_html=False, include_plotlyjs=False) if has_pca else ""
    
    # Build the HTML dashboard
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RNA-Seq Quality Control Dashboard</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: 'Inter', 'Segoe UI', system-ui, sans-serif;
            background: linear-gradient(135deg, {COLORS["background"]} 0%, #0f0f1a 100%);
            min-height: 100vh;
            color: {COLORS["text"]};
            padding: 20px;
        }}
        .dashboard-header {{
            text-align: center;
            padding: 30px 20px;
            margin-bottom: 30px;
            background: linear-gradient(135deg, {COLORS["surface"]} 0%, {COLORS["background"]} 100%);
            border-radius: 16px;
            border: 1px solid {COLORS["grid"]};
            box-shadow: 0 10px 40px rgba(0,0,0,0.3);
        }}
        .dashboard-header h1 {{
            font-size: 2.5rem;
            font-weight: 700;
            background: linear-gradient(90deg, {COLORS["primary"]}, {COLORS["secondary"]});
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 10px;
        }}
        .dashboard-header p {{
            color: #888;
            font-size: 1.1rem;
        }}
        .stats-summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .stat-card {{
            background: {COLORS["surface"]};
            border-radius: 12px;
            padding: 24px;
            text-align: center;
            border: 1px solid {COLORS["grid"]};
            transition: transform 0.3s, box-shadow 0.3s;
        }}
        .stat-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 10px 30px rgba(0, 166, 166, 0.2);
        }}
        .stat-value {{
            font-size: 2rem;
            font-weight: 700;
            color: {COLORS["primary"]};
        }}
        .stat-label {{
            color: #888;
            font-size: 0.9rem;
            margin-top: 5px;
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 20px;
            margin-bottom: 20px;
        }}
        @media (max-width: 1200px) {{
            .plot-grid {{
                grid-template-columns: 1fr;
            }}
        }}
        .plot-container {{
            background: {COLORS["surface"]};
            border-radius: 12px;
            padding: 15px;
            border: 1px solid {COLORS["grid"]};
        }}
        .plot-container.full-width {{
            grid-column: 1 / -1;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 0.85rem;
        }}
    </style>
</head>
<body>
    <div class="dashboard-header">
        <h1>RNA-Seq Quality Control Dashboard</h1>
        <p>Transcriptome Health Assessment | n = {len(dataset.sample_ids)} samples | {len(dataset.gene_ids):,} genes</p>
    </div>
    
    <div class="stats-summary">
        <div class="stat-card">
            <div class="stat-value">{len(dataset.sample_ids)}</div>
            <div class="stat-label">Total Samples</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{len(dataset.filtered_genes):,}</div>
            <div class="stat-label">Expressed Genes</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{qc_df["Library_Size"].mean()/1e6:.1f}M</div>
            <div class="stat-label">Mean Sequencing Depth</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{qc_df["Detected_Genes"].mean():,.0f}</div>
            <div class="stat-label">Mean Gene Detection</div>
        </div>
        <div class="stat-card">
            <div class="stat-value" style="color: {COLORS["warning"] if (qc_df["Library_Size"] < lib_threshold).sum() > 0 else COLORS["success"]}">{(qc_df["Library_Size"] < lib_threshold).sum()}</div>
            <div class="stat-label">QC-Failed Samples</div>
        </div>
    </div>
    
    <div class="plot-grid">
        <div class="plot-container">
            {lib_html}
        </div>
        <div class="plot-container">
            {genes_html}
        </div>
    </div>
    
    <div class="plot-grid">
        {"<div class='plot-container'>" + pca_html + "</div>" if has_pca else ""}
        <div class="plot-container">
            {mt_html}
        </div>
    </div>
    
    <div class="plot-container full-width" style="margin-top: 20px;">
        {summary_html}
    </div>
    
    <div class="footer">
        <p>Generated by Transcriptome Health Dashboard v2.0</p>
    </div>
</body>
</html>'''
    
    # Save dashboard
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logger.info(f"Dashboard saved to: {output_path}")
    
    return output_path


def save_individual_plots(
    dataset: 'RNASeqDataset',
    output_dir: str,
    lib_threshold: int = 20_000_000,
    mt_threshold: float = 20.0
) -> Dict[str, str]:
    """
    Save each QC plot as a separate HTML file.
    
    Args:
        dataset: RNASeqDataset instance
        output_dir: Directory to save plots
        lib_threshold: Library size threshold
        mt_threshold: MT percentage threshold
        
    Returns:
        Dictionary mapping plot names to file paths
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_files = {}
    
    qc_df = dataset.qc_metrics
    
    # Library sizes
    fig = plot_library_sizes_interactive(qc_df, threshold=lib_threshold)
    path = os.path.join(output_dir, "library_sizes.html")
    fig.write_html(path)
    saved_files["library_sizes"] = path
    
    # Detected genes
    fig = plot_detected_genes_interactive(qc_df)
    path = os.path.join(output_dir, "detected_genes.html")
    fig.write_html(path)
    saved_files["detected_genes"] = path
    
    # MT content
    fig = plot_mt_content_interactive(qc_df, threshold=mt_threshold)
    path = os.path.join(output_dir, "mt_content.html")
    fig.write_html(path)
    saved_files["mt_content"] = path
    
    # PCA if available
    if dataset.pca_results is not None:
        fig = plot_pca_interactive(
            dataset.pca_results["coordinates"],
            dataset.pca_results["explained_variance"],
            dataset.sample_ids
        )
        path = os.path.join(output_dir, "pca.html")
        fig.write_html(path)
        saved_files["pca"] = path
    
    return saved_files


# Legacy compatibility: Keep static plotting functions
def plot_library_sizes(qc_df: pd.DataFrame, output_dir: str) -> str:
    """Legacy static plot - redirects to interactive version."""
    fig = plot_library_sizes_interactive(qc_df)
    path = os.path.join(output_dir, "qc_library_sizes.html")
    fig.write_html(path)
    return path


    return path

def generate_responsive_figures(dataset: 'RNASeqDataset', output_dir: str):
    """Generate both HTML and static PNG figures for documentation."""
    os.makedirs(output_dir, exist_ok=True)
    
    if dataset.qc_metrics is None:
        dataset.calculate_qc_metrics()
    
    # Library Size
    fig = plot_library_sizes_interactive(dataset.qc_metrics)
    fig.write_html(os.path.join(output_dir, "library_size.html"))
    try:
        fig.write_image(os.path.join(output_dir, "library_size.png"), scale=2)
    except Exception as e:
        print(f"Could not save static image: {e}")

    # Detected Genes
    fig = plot_detected_genes_interactive(dataset.qc_metrics)
    fig.write_html(os.path.join(output_dir, "detected_genes.html"))
    try:
        fig.write_image(os.path.join(output_dir, "detected_genes.png"), scale=2)
    except Exception as e:
        print(f"Could not save static image: {e}")

    # MT Content - ensure function exists or use interactive one? 
    # The interactive one is plot_mt_content_interactive
    fig = plot_mt_content_interactive(dataset.qc_metrics)
    fig.write_html(os.path.join(output_dir, "mt_content.html"))
    try:
        fig.write_image(os.path.join(output_dir, "mt_content.png"), scale=2)
    except Exception as e:
        print(f"Could not save static image: {e}")

    # PCA
    if dataset.pca_results:
        fig = plot_pca_interactive(
            dataset.pca_results["coordinates"],
            dataset.pca_results["explained_variance"],
            dataset.sample_ids
        )
        fig.write_html(os.path.join(output_dir, "pca.html"))
        try:
            fig.write_image(os.path.join(output_dir, "pca.png"), scale=2)
        except Exception as e:
            print(f"Could not save static image: {e}")

if __name__ == "__main__":
    from dataset import RNASeqDataset
    # Mock data loading for visualization generation if needed, or load real data
    # Assuming we run this from root with data/TCGA...
    data_path = "data/TCGA_LIHC_Gene_Expression.csv"
    if os.path.exists(data_path):
        ds = RNASeqDataset(data_path)
        print("Loading data...")
        ds.load_data(data_path)
        print("Generating figures...")
        generate_responsive_figures(ds, "results")
        print("Done.")
    else:
        print("Data file not found.")