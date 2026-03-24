"""
Complete Sensitivity Analysis Script
Generates 9 types of charts for sensitivity analysis of hydrological models
Author: Dynamic Analysis
Date: 2025-11-05
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, pearsonr
import warnings
warnings.filterwarnings('ignore')

# Chart style configuration
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

class SensitivityAnalyzer:
    """
    Main class for sensitivity analysis with multiple chart types
    """

    def __init__(self, file_path, output_dir='sensitivity_results'):
        """
        Initialize the analyzer

        Args:
            file_path: Path to CSV file with parameters and simulations
            output_dir: Directory where charts will be saved
        """
        self.file_path = file_path
        self.output_dir = output_dir
        self.df = None
        self.parameters = []
        self.simulations = []
        self.n_params = 0
        self.n_sims = 0

        # Create output directory if it doesn't exist
        import os
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def load_data(self):
        """
        Loads data from CSV file dynamically
        Automatically identifies parameter and simulation columns
        """
        print("="*100)
        print("LOADING SENSITIVITY ANALYSIS DATA")
        print("="*100)

        # Load data
        self.df = pd.read_csv(self.file_path)
        print(f"✓ File loaded: {self.file_path}")
        print(f"✓ Dimensions: {self.df.shape[0]} rows × {self.df.shape[1]} columns")

        # Identify parameter columns (don't start with 'simulation')
        self.parameters = [col for col in self.df.columns if not col.startswith('simulation')]
        self.n_params = len(self.parameters)
        print(f"✓ Parameters identified: {self.n_params}")
        print(f"  Parameters: {', '.join(self.parameters)}")

        # Identify simulation columns
        self.simulations = [col for col in self.df.columns if col.startswith('simulation')]
        self.n_sims = len(self.simulations)
        print(f"✓ Simulations identified: {self.n_sims}")

        # Basic statistics
        print(f"\n📊 BASIC STATISTICS:")
        print(f"  Parameter ranges:")
        for param in self.parameters:
            print(f"    {param:20s}: [{self.df[param].min():10.4f}, {self.df[param].max():10.4f}]")

        print(f"\n  Simulation ranges:")
        sim_values = self.df[self.simulations].values.flatten()
        print(f"    Minimum: {sim_values.min():10.4f}")
        print(f"    Maximum: {sim_values.max():10.4f}")
        print(f"    Mean:    {sim_values.mean():10.4f}")
        print(f"    Median:  {np.median(sim_values):10.4f}")

        return self

    def chart_1_parameter_boxplots(self):
        """
        CHART 1: Boxplots of parameter distribution
        Shows the distribution of values for each parameter
        """
        print(f"\n{'='*100}")
        print("CHART 1: PARAMETER BOXPLOTS")
        print("="*100)

        n_params = len(self.parameters)
        n_cols = 4
        n_rows = int(np.ceil(n_params / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, n_rows*4))
        axes = axes.flatten() if n_params > 1 else [axes]

        for idx, param in enumerate(self.parameters):
            ax = axes[idx]

            # Boxplot with points
            bp = ax.boxplot([self.df[param]], vert=True, patch_artist=True,
                           boxprops=dict(facecolor='lightblue', alpha=0.7),
                           medianprops=dict(color='red', linewidth=2),
                           whiskerprops=dict(linewidth=1.5),
                           capprops=dict(linewidth=1.5))

            # Add individual points
            y = self.df[param]
            x = np.random.normal(1, 0.04, size=len(y))
            ax.scatter(x, y, alpha=0.3, s=20, c='navy')

            ax.set_title(f'{param}', fontsize=12, fontweight='bold')
            ax.set_ylabel('Value', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.set_xticks([])

            # Statistics
            stats_text = f'μ={y.mean():.2f}\nσ={y.std():.2f}'
            ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
                   verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                   fontsize=8)

        # Hide extra axes
        for idx in range(n_params, len(axes)):
            axes[idx].set_visible(False)

        plt.suptitle('MODEL PARAMETER DISTRIBUTION',
                    fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()

        output_file = f'{self.output_dir}/01_parameter_boxplots.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

    def chart_2_scatter_parameters_vs_simulations(self):
        """
        CHART 2: Scatter plots of parameters vs simulation results
        Shows the relationship between each parameter and mean of simulations
        """
        print(f"\n{'='*100}")
        print("CHART 2: SCATTER PLOTS PARAMETERS VS SIMULATIONS")
        print("="*100)

        # Calculate simulation statistics
        self.df['sim_mean'] = self.df[self.simulations].mean(axis=1)
        self.df['sim_std'] = self.df[self.simulations].std(axis=1)
        self.df['sim_cv'] = self.df['sim_std'] / self.df['sim_mean']  # Coefficient of variation

        n_params = len(self.parameters)
        n_cols = 4
        n_rows = int(np.ceil(n_params / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, n_rows*4))
        axes = axes.flatten() if n_params > 1 else [axes]

        for idx, param in enumerate(self.parameters):
            ax = axes[idx]

            x = self.df[param]
            y = self.df['sim_mean']

            # Scatter plot with colormap by density
            ax.scatter(x, y, alpha=0.6, s=30, c=self.df['sim_cv'],
                      cmap='viridis', edgecolors='black', linewidth=0.5)

            # Trend line
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            ax.plot(x, p(x), "r--", alpha=0.8, linewidth=2, label='Linear trend')

            # Calculate correlation
            corr_pearson, p_value = pearsonr(x, y)
            corr_spearman, _ = spearmanr(x, y)

            ax.set_xlabel(param, fontsize=10, fontweight='bold')
            ax.set_ylabel('Simulation Mean', fontsize=10)
            ax.set_title(f'{param}', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)

            # Correlation statistics
            stats_text = f'Pearson: {corr_pearson:.3f}\nSpearman: {corr_spearman:.3f}'
            color = 'green' if abs(corr_pearson) > 0.5 else 'orange' if abs(corr_pearson) > 0.3 else 'red'
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor=color, alpha=0.3),
                   fontsize=8)

        # Hide extra axes
        for idx in range(n_params, len(axes)):
            axes[idx].set_visible(False)

        plt.suptitle('PARAMETERS VS SIMULATION MEAN RELATIONSHIP\n(Color indicates coefficient of variation)',
                    fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()

        output_file = f'{self.output_dir}/02_scatter_parameters_vs_simulations.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

    def chart_3_correlation_heatmap(self):
        """
        CHART 3: Correlation heatmap between parameters and simulation statistics
        """
        print(f"\n{'='*100}")
        print("CHART 3: CORRELATION HEATMAP")
        print("="*100)

        # Prepare data for correlation
        simulation_stats = pd.DataFrame({
            'Sim_Mean': self.df[self.simulations].mean(axis=1),
            'Sim_Median': self.df[self.simulations].median(axis=1),
            'Sim_Std': self.df[self.simulations].std(axis=1),
            'Sim_Min': self.df[self.simulations].min(axis=1),
            'Sim_Max': self.df[self.simulations].max(axis=1),
            'Sim_Q25': self.df[self.simulations].quantile(0.25, axis=1),
            'Sim_Q75': self.df[self.simulations].quantile(0.75, axis=1)
        })

        # Combine parameters with statistics
        df_correlation = pd.concat([self.df[self.parameters], simulation_stats], axis=1)

        # Calculate correlation matrix
        corr_matrix = df_correlation.corr()

        # Create heatmap
        fig, ax = plt.subplots(figsize=(16, 14))

        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='RdBu_r',
                   center=0, vmin=-1, vmax=1,
                   square=True, linewidths=0.5,
                   cbar_kws={"shrink": 0.8, "label": "Pearson Correlation"},
                   ax=ax)

        ax.set_title('CORRELATION MATRIX: PARAMETERS AND SIMULATION STATISTICS',
                    fontsize=16, fontweight='bold', pad=20)

        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()

        output_file = f'{self.output_dir}/03_correlation_heatmap.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

    def chart_4_histograms(self):
        """
        CHART 4: Histograms of parameter and simulation distribution
        """
        print(f"\n{'='*100}")
        print("CHART 4: DISTRIBUTION HISTOGRAMS")
        print("="*100)

        # Figure 4a: Parameter histograms
        n_params = len(self.parameters)
        n_cols = 4
        n_rows = int(np.ceil(n_params / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, n_rows*4))
        axes = axes.flatten() if n_params > 1 else [axes]

        for idx, param in enumerate(self.parameters):
            ax = axes[idx]
            data = self.df[param]

            # Histogram
            n, bins, patches = ax.hist(data, bins=30, alpha=0.7, color='skyblue',
                                       edgecolor='black', density=True)

            # Density curve
            mu, std = data.mean(), data.std()
            xmin, xmax = ax.get_xlim()
            x = np.linspace(xmin, xmax, 100)
            p = stats.norm.pdf(x, mu, std)
            ax.plot(x, p, 'r-', linewidth=2, label='Normal Distribution')

            ax.set_xlabel(param, fontsize=10, fontweight='bold')
            ax.set_ylabel('Density', fontsize=10)
            ax.set_title(f'{param}', fontsize=12, fontweight='bold')
            ax.legend(fontsize=8, loc='upper right', framealpha=0.9)
            ax.grid(True, alpha=0.3)

            # Normality test (Shapiro-Wilk)
            if len(data) < 5000:  # Shapiro-Wilk doesn't work well with very large samples
                _, p_value = stats.shapiro(data)
                normal_text = 'Normal' if p_value > 0.05 else 'Not Normal'
                color = 'green' if p_value > 0.05 else 'red'
                ax.text(0.98, 0.98, f'{normal_text}\np={p_value:.4f}',
                       transform=ax.transAxes, verticalalignment='top',
                       horizontalalignment='right',
                       bbox=dict(boxstyle='round', facecolor=color, alpha=0.3),
                       fontsize=8)

        # Hide extra axes
        for idx in range(n_params, len(axes)):
            axes[idx].set_visible(False)

        plt.suptitle('PARAMETER DISTRIBUTION WITH NORMAL FIT',
                    fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()

        output_file = f'{self.output_dir}/04a_parameter_histograms.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

        # Figure 4b: Histogram of all simulations
        fig, ax = plt.subplots(figsize=(12, 6))

        sim_values = self.df[self.simulations].values.flatten()

        n, bins, patches = ax.hist(sim_values, bins=50, alpha=0.7, color='lightcoral',
                                   edgecolor='black', density=True)

        # Density curve
        mu, std = sim_values.mean(), sim_values.std()
        xmin, xmax = ax.get_xlim()
        x = np.linspace(xmin, xmax, 100)
        p = stats.norm.pdf(x, mu, std)
        ax.plot(x, p, 'b-', linewidth=2, label='Normal Distribution')

        ax.set_xlabel('Simulation Value', fontsize=12, fontweight='bold')
        ax.set_ylabel('Density', fontsize=12, fontweight='bold')
        ax.set_title(f'DISTRIBUTION OF ALL SIMULATIONS (n={self.n_sims})',
                    fontsize=14, fontweight='bold')
        ax.legend(fontsize=10, loc='upper right', framealpha=0.9)
        ax.grid(True, alpha=0.3)

        # Statistics
        stats_text = f'Mean: {mu:.2f}\nMedian: {np.median(sim_values):.2f}\nStd Dev: {std:.2f}'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
               fontsize=10)

        plt.tight_layout()

        output_file = f'{self.output_dir}/04b_simulation_histogram.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

    def chart_5_time_series(self):
        """
        CHART 5: Time series of simulations
        Shows the evolution of simulations over time
        """
        print(f"\n{'='*100}")
        print("CHART 5: SIMULATION TIME SERIES")
        print("="*100)

        # Transpose to have simulations as index and rows as time series
        sim_data = self.df[self.simulations].T

        fig, axes = plt.subplots(3, 1, figsize=(16, 12))

        # Chart 5a: All series (limited sample if too many)
        ax = axes[0]
        n_samples = min(100, len(self.df))  # Maximum 100 series for visualization
        sample_indices = np.random.choice(len(self.df), n_samples, replace=False)

        for idx in sample_indices:
            ax.plot(sim_data.index, sim_data[idx], alpha=0.3, linewidth=0.5)

        ax.set_xlabel('Simulation Number', fontsize=10)
        ax.set_ylabel('Value', fontsize=10)
        ax.set_title(f'ALL TIME SERIES (Sample of {n_samples} series)',
                    fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Format x-axis to avoid label overlap
        ax.xaxis.set_major_locator(plt.MaxNLocator(10))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        # Chart 5b: Statistics by simulation
        ax = axes[1]

        mean_vals = sim_data.mean(axis=1)
        std_vals = sim_data.std(axis=1)

        ax.plot(sim_data.index, mean_vals, 'b-', linewidth=2, label='Mean')
        ax.fill_between(sim_data.index,
                        mean_vals - std_vals,
                        mean_vals + std_vals,
                        alpha=0.3, color='blue', label='± 1 Std Dev')

        ax.plot(sim_data.index, sim_data.median(axis=1), 'r--', linewidth=2, label='Median')

        ax.set_xlabel('Simulation Number', fontsize=10)
        ax.set_ylabel('Value', fontsize=10)
        ax.set_title('STATISTICS BY SIMULATION (Mean ± Std Dev)',
                    fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Format x-axis to avoid label overlap
        ax.xaxis.set_major_locator(plt.MaxNLocator(10))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        # Chart 5c: Percentiles
        ax = axes[2]

        p10 = sim_data.quantile(0.10, axis=1)
        p25 = sim_data.quantile(0.25, axis=1)
        p50 = sim_data.quantile(0.50, axis=1)
        p75 = sim_data.quantile(0.75, axis=1)
        p90 = sim_data.quantile(0.90, axis=1)

        ax.fill_between(sim_data.index, p10, p90, alpha=0.2, color='gray', label='P10-P90')
        ax.fill_between(sim_data.index, p25, p75, alpha=0.4, color='blue', label='P25-P75')
        ax.plot(sim_data.index, p50, 'r-', linewidth=2, label='Median')

        ax.set_xlabel('Simulation Number', fontsize=10)
        ax.set_ylabel('Value', fontsize=10)
        ax.set_title('SIMULATION PERCENTILES', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Format x-axis to avoid label overlap
        ax.xaxis.set_major_locator(plt.MaxNLocator(10))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        plt.suptitle('TIME SERIES ANALYSIS OF SIMULATIONS',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()

        output_file = f'{self.output_dir}/05_time_series.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

    def chart_6_sensitivity_indices(self):
        """
        CHART 6: Sensitivity analysis with indices
        Calculates sensitivity indices based on correlation and regression
        """
        print(f"\n{'='*100}")
        print("CHART 6: SENSITIVITY INDICES")
        print("="*100)

        # Calculate different sensitivity metrics
        sensitivity = {
            'Parameter': [],
            'Pearson_Correlation': [],
            'Spearman_Correlation': [],
            'R2_Regression': [],
            'Sensitivity_Index': []
        }

        sim_mean = self.df[self.simulations].mean(axis=1)

        for param in self.parameters:
            x = self.df[param]

            # Pearson correlation
            corr_p, _ = pearsonr(x, sim_mean)

            # Spearman correlation
            corr_s, _ = spearmanr(x, sim_mean)

            # R² from linear regression
            from sklearn.linear_model import LinearRegression
            model = LinearRegression()
            model.fit(x.values.reshape(-1, 1), sim_mean)
            r2 = model.score(x.values.reshape(-1, 1), sim_mean)

            # Normalized sensitivity index (average of absolute metrics)
            index = (abs(corr_p) + abs(corr_s) + r2) / 3

            sensitivity['Parameter'].append(param)
            sensitivity['Pearson_Correlation'].append(corr_p)
            sensitivity['Spearman_Correlation'].append(corr_s)
            sensitivity['R2_Regression'].append(r2)
            sensitivity['Sensitivity_Index'].append(index)

        df_sens = pd.DataFrame(sensitivity)
        df_sens = df_sens.sort_values('Sensitivity_Index', ascending=False)

        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # Chart 6a: Sensitivity index
        ax = axes[0, 0]
        colors = plt.cm.RdYlGn(df_sens['Sensitivity_Index'] / df_sens['Sensitivity_Index'].max())
        bars = ax.barh(df_sens['Parameter'], df_sens['Sensitivity_Index'], color=colors)
        ax.set_xlabel('Normalized Sensitivity Index', fontsize=10, fontweight='bold')
        ax.set_title('PARAMETER SENSITIVITY RANKING', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')

        # Add values on bars
        for i, (param, val) in enumerate(zip(df_sens['Parameter'], df_sens['Sensitivity_Index'])):
            ax.text(val, i, f' {val:.3f}', va='center', fontsize=8)

        # Chart 6b: Pearson correlation
        ax = axes[0, 1]
        colors_p = ['green' if x > 0 else 'red' for x in df_sens['Pearson_Correlation']]
        ax.barh(df_sens['Parameter'], df_sens['Pearson_Correlation'], color=colors_p, alpha=0.7)
        ax.set_xlabel('Pearson Correlation', fontsize=10, fontweight='bold')
        ax.set_title('LINEAR CORRELATION (Pearson)', fontsize=12, fontweight='bold')
        ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
        ax.grid(True, alpha=0.3, axis='x')

        # Chart 6c: Spearman correlation
        ax = axes[1, 0]
        colors_s = ['green' if x > 0 else 'red' for x in df_sens['Spearman_Correlation']]
        ax.barh(df_sens['Parameter'], df_sens['Spearman_Correlation'], color=colors_s, alpha=0.7)
        ax.set_xlabel('Spearman Correlation', fontsize=10, fontweight='bold')
        ax.set_title('NON-LINEAR CORRELATION (Spearman)', fontsize=12, fontweight='bold')
        ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
        ax.grid(True, alpha=0.3, axis='x')

        # Chart 6d: R² regression
        ax = axes[1, 1]
        colors_r2 = plt.cm.YlOrRd(df_sens['R2_Regression'])
        ax.barh(df_sens['Parameter'], df_sens['R2_Regression'], color=colors_r2)
        ax.set_xlabel('R² Linear Regression', fontsize=10, fontweight='bold')
        ax.set_title('LINEAR FIT GOODNESS', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')

        plt.suptitle('SENSITIVITY ANALYSIS: MULTIPLE INDICES',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()

        output_file = f'{self.output_dir}/06_sensitivity_indices.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")

        # Save sensitivity table
        csv_file = f'{self.output_dir}/06_sensitivity_table.csv'
        df_sens.to_csv(csv_file, index=False)
        print(f"✓ Table saved: {csv_file}")

        plt.close()

    def chart_7_violin_plots(self):
        """
        CHART 7: Violin plots of parameters
        Combines boxplot with probability density
        """
        print(f"\n{'='*100}")
        print("CHART 7: VIOLIN PLOTS")
        print("="*100)

        # Prepare data in long format for seaborn
        df_long = self.df[self.parameters].melt(var_name='Parameter', value_name='Value')

        fig, ax = plt.subplots(figsize=(16, 8))

        # Violin plot
        sns.violinplot(data=df_long, x='Parameter', y='Value', ax=ax,
                      palette='Set2', inner='box', cut=0)

        ax.set_xlabel('Parameters', fontsize=12, fontweight='bold')
        ax.set_ylabel('Value', fontsize=12, fontweight='bold')
        ax.set_title('PARAMETER DISTRIBUTION (Violin Plots)',
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        output_file = f'{self.output_dir}/07_violin_plots.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

    def chart_8_parallel_coordinates(self):
        """
        CHART 8: Parallel coordinates plot
        Multidimensional visualization of parameters
        """
        print(f"\n{'='*100}")
        print("CHART 8: PARALLEL COORDINATES")
        print("="*100)

        from pandas.plotting import parallel_coordinates

        # Normalize parameters for visualization
        df_norm = self.df[self.parameters].copy()
        for col in self.parameters:
            min_val = df_norm[col].min()
            max_val = df_norm[col].max()
            df_norm[col] = (df_norm[col] - min_val) / (max_val - min_val)

        # Add category based on simulation mean
        sim_mean = self.df[self.simulations].mean(axis=1)
        df_norm['Category'] = pd.cut(sim_mean, bins=5, labels=['Very Low', 'Low', 'Medium', 'High', 'Very High'])

        # Take sample if too many rows
        n_sample = min(500, len(df_norm))
        df_sample = df_norm.sample(n=n_sample, random_state=42)

        fig, ax = plt.subplots(figsize=(18, 8))

        parallel_coordinates(df_sample, 'Category', ax=ax,
                           colormap='viridis', alpha=0.5, linewidth=0.8)

        ax.set_xlabel('Parameters', fontsize=12, fontweight='bold')
        ax.set_ylabel('Normalized Value [0-1]', fontsize=12, fontweight='bold')
        ax.set_title(f'PARALLEL COORDINATES (Sample of {n_sample} cases)\nColored by simulation result',
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        ax.legend(title='Simulation Result', loc='upper left', fontsize=10)

        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()

        output_file = f'{self.output_dir}/08_parallel_coordinates.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

    def chart_9_tornado_plot(self):
        """
        CHART 9: Tornado plot
        Shows the impact of each parameter on the range of results
        """
        print(f"\n{'='*100}")
        print("CHART 9: TORNADO PLOT")
        print("="*100)

        # Calculate result range for each parameter
        tornado_data = {
            'Parameter': [],
            'Range_Min': [],
            'Range_Max': [],
            'Amplitude': []
        }

        sim_mean = self.df[self.simulations].mean(axis=1)
        baseline = sim_mean.mean()

        for param in self.parameters:
            # Sort by parameter and take extremes
            df_sorted = self.df.sort_values(param)

            # Take bottom and top 10%
            n_sample = max(1, int(len(df_sorted) * 0.1))

            low_result = df_sorted.head(n_sample)[self.simulations].mean().mean()
            high_result = df_sorted.tail(n_sample)[self.simulations].mean().mean()

            amplitude = abs(high_result - low_result)

            tornado_data['Parameter'].append(param)
            tornado_data['Range_Min'].append(min(low_result, high_result) - baseline)
            tornado_data['Range_Max'].append(max(low_result, high_result) - baseline)
            tornado_data['Amplitude'].append(amplitude)

        df_tornado = pd.DataFrame(tornado_data)
        df_tornado = df_tornado.sort_values('Amplitude', ascending=True)

        # Create chart
        fig, ax = plt.subplots(figsize=(12, 10))

        y_pos = np.arange(len(df_tornado))

        # Bars to the left (negative) and right (positive)
        left_bars = ax.barh(y_pos, df_tornado['Range_Min'],
                           color='coral', alpha=0.7, label='Decrease')
        right_bars = ax.barh(y_pos, df_tornado['Range_Max'],
                            color='skyblue', alpha=0.7, label='Increase')

        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_tornado['Parameter'])
        ax.set_xlabel('Change in Result (relative to baseline)', fontsize=12, fontweight='bold')
        ax.set_title('TORNADO PLOT: PARAMETER IMPACT ON RESULTS\n(Sorted by impact amplitude)',
                    fontsize=14, fontweight='bold')
        ax.axvline(x=0, color='black', linestyle='-', linewidth=1.5)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='x')

        # Add amplitude values
        for i, (param, amp) in enumerate(zip(df_tornado['Parameter'], df_tornado['Amplitude'])):
            ax.text(0, i, f' {amp:.1f} ', va='center', ha='center',
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5),
                   fontsize=8)

        plt.tight_layout()

        output_file = f'{self.output_dir}/09_tornado_plot.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Chart saved: {output_file}")
        plt.close()

    def generate_summary_report(self):
        """
        Generates a summary report in text format
        """
        print(f"\n{'='*100}")
        print("GENERATING SUMMARY REPORT")
        print("="*100)

        report_file = f'{self.output_dir}/00_summary_report.txt'

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*100 + "\n")
            f.write("SENSITIVITY ANALYSIS REPORT\n")
            f.write("="*100 + "\n\n")

            f.write(f"File analyzed: {self.file_path}\n")
            f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            f.write("DATA SUMMARY:\n")
            f.write("-"*100 + "\n")
            f.write(f"  • Number of parameter sets: {len(self.df)}\n")
            f.write(f"  • Number of parameters: {self.n_params}\n")
            f.write(f"  • Number of simulations: {self.n_sims}\n\n")

            f.write("MODEL PARAMETERS:\n")
            f.write("-"*100 + "\n")
            for param in self.parameters:
                f.write(f"  • {param:20s}: ")
                f.write(f"[{self.df[param].min():10.4f}, {self.df[param].max():10.4f}] ")
                f.write(f"(μ={self.df[param].mean():10.4f}, σ={self.df[param].std():10.4f})\n")

            f.write("\n")
            f.write("SIMULATION STATISTICS:\n")
            f.write("-"*100 + "\n")
            sim_values = self.df[self.simulations].values.flatten()
            f.write(f"  • Minimum:  {sim_values.min():15.4f}\n")
            f.write(f"  • Maximum:  {sim_values.max():15.4f}\n")
            f.write(f"  • Mean:     {sim_values.mean():15.4f}\n")
            f.write(f"  • Median:   {np.median(sim_values):15.4f}\n")
            f.write(f"  • Std Dev:  {sim_values.std():15.4f}\n")
            f.write(f"  • Q25:      {np.percentile(sim_values, 25):15.4f}\n")
            f.write(f"  • Q75:      {np.percentile(sim_values, 75):15.4f}\n")

            f.write("\n")
            f.write("CHARTS GENERATED:\n")
            f.write("-"*100 + "\n")
            f.write("  01. Parameter distribution boxplots\n")
            f.write("  02. Scatter plots: Parameters vs Simulations\n")
            f.write("  03. Correlation heatmap\n")
            f.write("  04. Distribution histograms\n")
            f.write("  05. Simulation time series\n")
            f.write("  06. Sensitivity indices\n")
            f.write("  07. Violin plots\n")
            f.write("  08. Parallel coordinates\n")
            f.write("  09. Tornado plot\n")

            f.write("\n" + "="*100 + "\n")
            f.write("END OF REPORT\n")
            f.write("="*100 + "\n")

        print(f"✓ Report saved: {report_file}")

    def run_complete_analysis(self):
        """
        Runs the complete analysis generating all charts
        """
        print("\n")
        print("╔" + "="*98 + "╗")
        print("║" + " "*32 + "COMPLETE SENSITIVITY ANALYSIS" + " "*37 + "║")
        print("╚" + "="*98 + "╝")
        print()

        # Load data
        self.load_data()

        # Generate all charts
        self.chart_1_parameter_boxplots()
        self.chart_2_scatter_parameters_vs_simulations()
        self.chart_3_correlation_heatmap()
        self.chart_4_histograms()
        self.chart_5_time_series()
        self.chart_6_sensitivity_indices()
        self.chart_7_violin_plots()
        self.chart_8_parallel_coordinates()
        self.chart_9_tornado_plot()

        # Generate report
        self.generate_summary_report()

        print(f"\n{'='*100}")
        print("✅ ANALYSIS COMPLETED SUCCESSFULLY")
        print("="*100)
        print(f"📁 All files saved in: {self.output_dir}/")
        print()


# =============================================================================
# MAIN FUNCTION TO RUN THE ANALYSIS
# =============================================================================
def main():
    """
    Main function to run the analysis
    """
    import sys

    # Configuration
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        # Use default file
        input_file = 'Params_Evaluation_Hy.csv'

    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    else:
        output_dir = 'sensitivity_results'

    print(f"\n📊 Starting sensitivity analysis...")
    print(f"📂 Input file: {input_file}")
    print(f"📁 Output directory: {output_dir}")

    # Create analyzer and run
    analyzer = SensitivityAnalyzer(input_file, output_dir)
    analyzer.run_complete_analysis()


if __name__ == "__main__":
    main()
