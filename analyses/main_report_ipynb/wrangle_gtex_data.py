

class gtexData:
    def __init__(
        self,
        df_metadata = None,
        df_tpm = None,
        df_counts = None
    ):
        self.df_metadata = df_metadata
        self.df_tpm = df_tpm.T
        self.df_counts = df_counts.T


    def data_shapes(self):
        """Print the shape of each array/dataframe in the object."""
        for attr in ['df_metadata', 'df_tpm', 'df_counts']:
            arr = getattr(self, attr)
            if arr is not None:
                print(f"{attr}: {arr.shape}")
            else:
                print(f"{attr}: None")

    def list_common_samples(self):
        """Return a list of row indices present in all non-None dataframes."""
        dfs = [self.df_metadata, self.df_tpm, self.df_counts]
        index_sets = [set(df.index) for df in dfs if df is not None]
        if not index_sets:
            return []
        common_indices = set.intersection(*index_sets)
        print(f"Number of common cell lines: {len(common_indices)}")
        return list(common_indices)
    
    def subset_to_common_samples(self):
        """Return a new gtexData object with dataframes subsetted to common samples."""
        common_indices = self.list_common_samples()
        def subset(df):
            if df is None:
                return None
            if hasattr(df, "loc"):  # pandas DataFrame
                return df.loc[common_indices]
            elif hasattr(df, "obs_names"):  # AnnData
                return df[common_indices, :]
            else:
                return df
        return gtexData(
            df_metadata=subset(self.df_metadata),
            df_tpm=subset(self.df_tpm),
            df_counts=subset(self.df_counts)
            
        )

    def subset_to_random_samples(self, n_samples=100, random_state=42):
        """
        Return a new gtexData object with dataframes subsetted to randomly selected samples.

        Parameters:
        -----------
        n_cells : int, default=100
            Number of cells to randomly select
        random_state : int, default=42
            Random seed for reproducibility
            
        Returns:
        --------
        gtexData object with randomly selected cells
        """
        import random
        
        # Get all available cell lines from all datasets
        all_indices = set()
        dfs = [self.df_metadata, self.df_tpm, self.df_counts]
        
        for df in dfs:
            if df is not None:
                if hasattr(df, "index"):  # pandas DataFrame
                    all_indices.update(df.index)
                elif hasattr(df, "obs_names"):  # AnnData
                    all_indices.update(df.obs_names)
        
        # Convert to list and randomly sample
        all_indices = list(all_indices)
        random.seed(random_state)
        
        # Select minimum of requested cells or available cells
        n_to_select = min(n_samples, len(all_indices))
        random_indices = random.sample(all_indices, n_to_select)

        print(f"Randomly selected {n_to_select} samples from {len(all_indices)} available samples")

        def subset(df):
            if df is None:
                return None
            if hasattr(df, "loc"):  # pandas DataFrame
                # Only select indices that exist in this dataframe
                available_indices = [idx for idx in random_indices if idx in df.index]
                return df.loc[available_indices] if available_indices else None
            elif hasattr(df, "obs_names"):  # AnnData
                # Only select indices that exist in this AnnData object
                available_indices = [idx for idx in random_indices if idx in df.obs_names]
                return df[available_indices, :] if available_indices else None
            else:
                return df
                
        return gtexData(
            df_metadata=subset(self.df_metadata),
            df_tpm=subset(self.df_tpm),
            df_counts=subset(self.df_counts)
        )
    
    def list_common_features(self):
        """Return a list of column names (features) present in all non-None objects (DataFrame or AnnData)."""
        dfs = [self.df_tpm, self.df_counts]
        feature_sets = []
        for df in dfs:
            if df is not None:
                if hasattr(df, "columns"):  # pandas DataFrame
                    feature_sets.append(set(df.columns))
                elif hasattr(df, "var_names"):  # AnnData
                    feature_sets.append(set(df.var_names))
        if not feature_sets:
            return []
        common_features = set.intersection(*feature_sets)
        print(f"Number of common features: {len(common_features)}")
        return list(common_features)
    
    def subset_to_common_features(self):
        """Return a new gtexData object with dataframes subsetted to common features (columns)."""
        common_features = self.list_common_features()
        def subset(df):
            if df is None:
                return None
            if hasattr(df, "loc"):  # pandas DataFrame
                return df.loc[:, common_features]
            elif hasattr(df, "var_names"):  # AnnData
                return df[:, common_features]
            else:
                return df
        return gtexData(
            df_metadata=self.df_metadata,
            df_tpm=subset(self.df_tpm),
            df_counts=subset(self.df_counts),
        )
    
    def feature_presence(self, feature_list):
        """
        For each feature in feature_list, return a dict indicating in which dataframes the feature is present.
        Returns: dict of {feature: [list of dataframe names where present]}
        """
        dfs = {
            'df_tpm': self.df_tpm,
            'df_counts': self.df_counts
        }
        result = {}
        for feature in feature_list:
            present_in = []
            for name, df in dfs.items():
                if df is not None:
                    if hasattr(df, "columns") and feature in df.columns:
                        present_in.append(name)
                    elif hasattr(df, "var_names") and feature in df.var_names:
                        present_in.append(name)
            result[feature] = present_in
        return result
    
    def save_data(self, filepath):
        """Save the gtexData object using pickle"""
        import pickle
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
    
    @staticmethod
    def load_data(filepath):
        """Load a gtexData object from file"""
        import pickle
        with open(filepath, 'rb') as f:
            return pickle.load(f)

## Done creating depmap data object                                             ##
##################################################################################
