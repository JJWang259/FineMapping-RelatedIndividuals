/*
 * LD Adjuster - Adjust LD Matrix with GRM and Covariates
 * 
 * Author: Jicai Jiang
 * Version: 1.1.1
 * 
 * This program adjusts linkage disequilibrium matrix using genomic 
 * relationship matrix and optional covariates to account for population 
 * structure, relatedness, and confounding factors in fine-mapping analysis.
 */

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Eigen;

// Use single precision as requested
typedef Matrix<float, Dynamic, Dynamic> MatrixXf_dyn;
typedef Matrix<float, Dynamic, 1> VectorXf_dyn;

class CSVReader {
public:
    vector<string> individual_ids;
    vector<string> column_names;
    unordered_map<string, vector<float>> data; // column_name -> values
    
    bool read(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot open CSV file: " << filename << endl;
            return false;
        }
        
        string line;
        
        // Read header
        if (!getline(file, line)) {
            cerr << "Error: Empty CSV file: " << filename << endl;
            return false;
        }
        
        istringstream header_iss(line);
        string token;
        vector<string> header;
        
        while (getline(header_iss, token, ',')) {
            // Trim whitespace
            token.erase(0, token.find_first_not_of(" \t\r\n"));
            token.erase(token.find_last_not_of(" \t\r\n") + 1);
            header.push_back(token);
        }
        
        if (header.empty()) {
            cerr << "Error: No columns found in CSV header" << endl;
            return false;
        }
        
        // Store column names (excluding first column which is IID)
        for (size_t i = 1; i < header.size(); ++i) {
            column_names.push_back(header[i]);
            data[header[i]] = vector<float>();
        }
        
        cout << "Found " << column_names.size() << " data columns in CSV file" << endl;
        
        // Read data rows
        while (getline(file, line)) {
            istringstream iss(line);
            vector<string> fields;
            
            while (getline(iss, token, ',')) {
                // Trim whitespace
                token.erase(0, token.find_first_not_of(" \t\r\n"));
                token.erase(token.find_last_not_of(" \t\r\n") + 1);
                fields.push_back(token);
            }
            
            if (fields.size() != header.size()) {
                cerr << "Warning: Skipping malformed CSV line (expected " 
                     << header.size() << " fields, got " << fields.size() << ")" << endl;
                continue;
            }
            
            // Store individual ID
            individual_ids.push_back(fields[0]);
            
            // Store data values
            bool has_missing = false;
            for (size_t i = 1; i < fields.size(); ++i) {
                try {
                    float val = stof(fields[i]);
                    data[header[i]].push_back(val);
                } catch (const invalid_argument&) {
                    has_missing = true;
                    break;
                }
            }
            
            if (has_missing) {
                // Remove this individual due to missing data
                individual_ids.pop_back();
                for (size_t i = 1; i < fields.size(); ++i) {
                    if (!data[header[i]].empty()) {
                        data[header[i]].pop_back();
                    }
                }
            }
        }
        
        file.close();
        cout << "Read " << individual_ids.size() << " individuals from CSV file" << endl;
        return true;
    }
    
    bool hasColumn(const string& col_name) const {
        return data.find(col_name) != data.end();
    }
    
    vector<float> getColumn(const string& col_name) const {
        auto it = data.find(col_name);
        if (it != data.end()) {
            return it->second;
        }
        return vector<float>();
    }
};

class PLINKRawReader {
public:
    vector<string> individual_ids;
    vector<string> snp_ids;
    MatrixXf_dyn genotype_matrix;
    
    bool read(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot open PLINK raw file: " << filename << endl;
            return false;
        }
        
        string line;
        vector<vector<float>> data;
        
        // Read header to get SNP IDs
        if (getline(file, line)) {
            istringstream iss(line);
            string token;
            vector<string> header;
            
            while (iss >> token) {
                header.push_back(token);
            }
            
            // SNP IDs start from column 6 (0-indexed)
            for (size_t i = 6; i < header.size(); ++i) {
                snp_ids.push_back(header[i]);
            }
            
            cout << "Found " << snp_ids.size() << " SNPs in raw file" << endl;
        }
        
        // Read individual data
        while (getline(file, line)) {
            istringstream iss(line);
            string token;
            vector<string> fields;
            
            while (iss >> token) {
                fields.push_back(token);
            }
            
            if (fields.size() < 6 + snp_ids.size()) {
                cerr << "Warning: Skipping malformed line" << endl;
                continue;
            }
            
            // Store individual ID (column 1, 0-indexed)
            individual_ids.push_back(fields[1]);
            
            // Store genotype data (starting from column 6)
            vector<float> genotypes;
            for (size_t i = 6; i < fields.size(); ++i) {
                try {
                    float val = stof(fields[i]);
                    genotypes.push_back(val);
                } catch (const invalid_argument&) {
                    genotypes.push_back(NAN); // Handle missing data
                }
            }
            data.push_back(genotypes);
        }
        
        file.close();
        
        // Convert to Eigen matrix (individuals x SNPs)
        int n_individuals = data.size();
        int n_snps = snp_ids.size();
        
        genotype_matrix = MatrixXf_dyn(n_individuals, n_snps);
        for (int i = 0; i < n_individuals; ++i) {
            for (int j = 0; j < n_snps; ++j) {
                genotype_matrix(i, j) = data[i][j];
            }
        }
        
        cout << "Read " << n_individuals << " individuals and " << n_snps << " SNPs" << endl;
        return true;
    }
};

class GRMReader {
public:
    vector<string> individual_ids;
    MatrixXf_dyn grm_matrix;
    
    bool read(const string& prefix) {
        string iid_file = prefix + ".grm.iid";
        string bin_file = prefix + ".grm.bin";
        
        // Read IID file
        ifstream iid_stream(iid_file);
        if (!iid_stream.is_open()) {
            cerr << "Error: Cannot open IID file: " << iid_file << endl;
            return false;
        }
        
        string line;
        while (getline(iid_stream, line)) {
            istringstream iss(line);
            string iid;
            if (iss >> iid) {
                individual_ids.push_back(iid);
            }
        }
        iid_stream.close();
        
        // Read binary GRM file
        ifstream bin_stream(bin_file, ios::binary);
        if (!bin_stream.is_open()) {
            cerr << "Error: Cannot open binary GRM file: " << bin_file << endl;
            return false;
        }
        
        // Read number of individuals
        int32_t np;
        bin_stream.read(reinterpret_cast<char*>(&np), sizeof(int32_t));
        
        if (np != static_cast<int32_t>(individual_ids.size())) {
            cerr << "Error: Mismatch of sample size between .grm.iid and .grm.bin" << endl;
            return false;
        }
        
        // Read number of markers
        float nm;
        bin_stream.read(reinterpret_cast<char*>(&nm), sizeof(float));
        
        // Construct matrix and read column by column directly
        grm_matrix = MatrixXf_dyn::Zero(np, np);
        
        for (int j = 0; j < np; ++j) {  // For each column
            // Read elements from diagonal down to bottom of column j
            int n_elements_in_col = np - j;
            vector<float> col_data(n_elements_in_col);
            bin_stream.read(reinterpret_cast<char*>(col_data.data()), n_elements_in_col * sizeof(float));
            
            // Fill the column from diagonal down, dividing by nm
            for (int i = 0; i < n_elements_in_col; ++i) {
                grm_matrix(j + i, j) = col_data[i] / nm;
            }
        }
        bin_stream.close();

        // Fill upper triangle using Eigen's efficient triangular view
        grm_matrix.triangularView<StrictlyUpper>() = grm_matrix.adjoint();
        
        cout << "Read GRM with " << np << " individuals, constructed from " << nm << " markers" << endl;

        return true;
    }
};

class LDAdjuster {
private:
    MatrixXf_dyn raw_data;
    MatrixXf_dyn covariate_data;
    const MatrixXf_dyn* original_grm;  // Pointer to original GRM (no copy)
    vector<int> grm_subset_indices;    // Indices of matched individuals in original GRM
    vector<string> snp_ids;
    vector<string> matched_individuals;
    vector<float> error_weights;       // Diagonal elements for R matrix
    bool use_covariates;
    float h2;
    
public:
    LDAdjuster(float heritability) : use_covariates(false), h2(heritability) {}
    
    bool loadData(const PLINKRawReader& raw_reader, const GRMReader& grm_reader, 
                  const CSVReader* weight_reader = nullptr, const string& weight_col_name = "",
                  const CSVReader* covariate_reader = nullptr, const vector<string>& covariate_col_names = {}) {
        
        // Create mapping from individual IDs to indices
        unordered_map<string, int> raw_id_to_idx;
        unordered_map<string, int> grm_id_to_idx;
        unordered_map<string, int> weight_id_to_idx;
        unordered_map<string, int> covariate_id_to_idx;
        
        for (size_t i = 0; i < raw_reader.individual_ids.size(); ++i) {
            raw_id_to_idx[raw_reader.individual_ids[i]] = i;
        }
        
        for (size_t i = 0; i < grm_reader.individual_ids.size(); ++i) {
            grm_id_to_idx[grm_reader.individual_ids[i]] = i;
        }
        
        if (weight_reader) {
            for (size_t i = 0; i < weight_reader->individual_ids.size(); ++i) {
                weight_id_to_idx[weight_reader->individual_ids[i]] = i;
            }
        }
        
        if (covariate_reader) {
            for (size_t i = 0; i < covariate_reader->individual_ids.size(); ++i) {
                covariate_id_to_idx[covariate_reader->individual_ids[i]] = i;
            }
        }
        
        // Find intersection of individuals across all required files
        vector<int> raw_indices, grm_indices, weight_indices, covariate_indices;
        
        for (const auto& raw_id : raw_reader.individual_ids) {
            auto grm_it = grm_id_to_idx.find(raw_id);
            if (grm_it == grm_id_to_idx.end()) continue;
            
            // Check weight file if provided
            if (weight_reader) {
                auto weight_it = weight_id_to_idx.find(raw_id);
                if (weight_it == weight_id_to_idx.end()) continue;
                weight_indices.push_back(weight_it->second);
            }
            
            // Check covariate file if provided
            if (covariate_reader) {
                auto covariate_it = covariate_id_to_idx.find(raw_id);
                if (covariate_it == covariate_id_to_idx.end()) continue;
                covariate_indices.push_back(covariate_it->second);
            }
            
            // If we reach here, individual is in all required files
            raw_indices.push_back(raw_id_to_idx[raw_id]);
            grm_indices.push_back(grm_it->second);
            matched_individuals.push_back(raw_id);
        }
        
        cout << "Found " << matched_individuals.size() << " matched individuals across all files" << endl;
        
        if (matched_individuals.size() == 0) {
            cerr << "Error: No matching individuals found across all required files" << endl;
            return false;
        }
        
        // Extract raw data subset
        int n_matched = matched_individuals.size();
        int n_snps = raw_reader.snp_ids.size();
        
        raw_data = MatrixXf_dyn(n_matched, n_snps);
        for (int i = 0; i < n_matched; ++i) {
            raw_data.row(i) = raw_reader.genotype_matrix.row(raw_indices[i]);
        }
        
        // Extract error weights if provided
        if (weight_reader) {
            vector<float> weight_col = weight_reader->getColumn(weight_col_name);
            error_weights.resize(n_matched);
            for (int i = 0; i < n_matched; ++i) {
                error_weights[i] = weight_col[weight_indices[i]];
            }
            cout << "Loaded error weights for " << n_matched << " individuals" << endl;
        } else {
            // Default to 1.0 for all individuals
            error_weights.resize(n_matched, 1.0f);
        }
        
        // Extract covariate data if provided
        if (covariate_reader && !covariate_col_names.empty()) {
            use_covariates = true;
            int n_covariates = covariate_col_names.size();
            covariate_data = MatrixXf_dyn(n_matched, n_covariates);
            
            for (int j = 0; j < n_covariates; ++j) {
                vector<float> covariate_col = covariate_reader->getColumn(covariate_col_names[j]);
                for (int i = 0; i < n_matched; ++i) {
                    covariate_data(i, j) = covariate_col[covariate_indices[i]];
                }
            }
            
            // Check if covariate matrix is full column rank
            MatrixXf_dyn Q_test = covariate_data;
            ColPivHouseholderQR<MatrixXf_dyn> qr(Q_test);
            if (qr.rank() < n_covariates) {
                cerr << "Error: Covariate matrix is not full column rank (rank = " 
                     << qr.rank() << ", columns = " << n_covariates << ")" << endl;
                return false;
            }
            
            cout << "Loaded " << n_covariates << " covariates for " << n_matched 
                 << " individuals (full rank confirmed)" << endl;
        }
        
        // Store GRM indices instead of copying the matrix
        grm_subset_indices = grm_indices;
        
        // Store reference to original GRM (no copy)
        original_grm = &grm_reader.grm_matrix;
        
        // Store SNP IDs
        snp_ids = raw_reader.snp_ids;
        
        // Filter out SNPs with no variation
        filterInvariantSNPs();
        
        return true;
    }
    
private:
    void filterInvariantSNPs() {
        vector<int> valid_snp_indices;
        vector<string> valid_snp_ids;
        
        for (int j = 0; j < raw_data.cols(); ++j) {
            // Calculate variance for this SNP
            VectorXf_dyn snp_col = raw_data.col(j);
            
            // Remove NaN values for variance calculation
            vector<float> valid_values;
            for (int i = 0; i < snp_col.size(); ++i) {
                if (!isnan(snp_col(i))) {
                    valid_values.push_back(snp_col(i));
                }
            }
            
            if (valid_values.size() < 2) continue;
            
            // Calculate variance
            float mean = 0.0f;
            for (float val : valid_values) mean += val;
            mean /= valid_values.size();
            
            float variance = 0.0f;
            for (float val : valid_values) {
                variance += (val - mean) * (val - mean);
            }
            variance /= (valid_values.size() - 1);
            
            if (variance > 1e-6f) { // Keep SNPs with non-zero variance
                valid_snp_indices.push_back(j);
                valid_snp_ids.push_back(snp_ids[j]);
            }
        }
        
        cout << "Filtered " << (snp_ids.size() - valid_snp_ids.size()) 
             << " invariant SNPs, keeping " << valid_snp_ids.size() << " SNPs" << endl;
        
        // Create new matrix with only valid SNPs
        MatrixXf_dyn filtered_data(raw_data.rows(), valid_snp_indices.size());
        for (size_t i = 0; i < valid_snp_indices.size(); ++i) {
            filtered_data.col(i) = raw_data.col(valid_snp_indices[i]);
        }
        
        raw_data = filtered_data;
        snp_ids = valid_snp_ids;
    }
    
public:
    bool computeLDAdjustment(const string& output_prefix) {
        auto start_time = chrono::high_resolution_clock::now();
        
        cout << "Starting LD adjustment computation..." << endl;
        cout << "Matrix dimensions: " << raw_data.rows() << " individuals x " 
             << raw_data.cols() << " SNPs" << endl;
        
        if (use_covariates) {
            cout << "Using " << covariate_data.cols() << " covariates" << endl;
        }
        
        // Step 1: Extract GRM subset and modify diagonal with error weights
        int n_matched = grm_subset_indices.size();
        MatrixXf_dyn grm_subset(n_matched, n_matched);
        
        float eta = h2 / (1.0f - h2);
        
        // Extract subset and apply heritability adjustment + error weights
        for (int i = 0; i < n_matched; ++i) {
            for (int j = 0; j < n_matched; ++j) {
                grm_subset(i, j) = (*original_grm)(grm_subset_indices[i], grm_subset_indices[j]);
            }
        }
        grm_subset *= eta;
        grm_subset.diagonal() += Map<VectorXf_dyn>(error_weights.data(), error_weights.size());
        
        cout << "Extracted GRM subset and applied heritability adjustment (eta = " << eta << ")" << endl;
        
        // Step 2: Standardize raw data
        for (int j = 0; j < raw_data.cols(); ++j) {
            VectorXf_dyn col = raw_data.col(j);
            
            // Handle missing values - replace NaN with column mean
            vector<float> valid_values;
            for (int i = 0; i < col.size(); ++i) {
                if (!isnan(col(i))) {
                    valid_values.push_back(col(i));
                }
            }
            
            if (valid_values.empty()) continue;
            
            float mean = 0.0f;
            for (float val : valid_values) mean += val;
            mean /= valid_values.size();
            
            // Replace NaN with mean and center
            for (int i = 0; i < col.size(); ++i) {
                if (isnan(raw_data(i, j))) {
                    raw_data(i, j) = mean;
                }
                raw_data(i, j) -= mean;
            }
            
            // Scale (compute std dev)
            float variance = raw_data.col(j).array().square().sum() / (raw_data.rows() - 1);
            float std_dev = sqrt(variance);
            
            if (std_dev > 1e-6f) {
                raw_data.col(j) /= std_dev;
            }
        }
        
        cout << "Standardized genotype matrix" << endl;
        
        // Step 3: Cholesky decomposition (in-place)
        LLT<MatrixXf_dyn> llt(grm_subset);
        if (llt.info() != Success) {
            cerr << "Error: Cholesky decomposition failed. Matrix may not be positive definite." << endl;
            return false;
        }
        
        cout << "Computed Cholesky decomposition" << endl;
        
        // Step 4: Solve L * X_star = X
        MatrixXf_dyn X_star = llt.matrixL().solve(raw_data);
        
        cout << "Solved linear system for genotypes" << endl;
        
        // Step 5: Handle covariates if present
        MatrixXf_dyn Q_star;
        if (use_covariates) {
            // Solve L * Q_star = Q
            Q_star = llt.matrixL().solve(covariate_data);
            cout << "Solved linear system for covariates" << endl;
        }
        
        // Step 6: Compute correlation matrix
        MatrixXf_dyn R_star;
        if (use_covariates) {
            // R* = (X*)^T X* - (X*)^T Q* ((Q*)^T Q*)^(-1) (Q*)^T X*
            MatrixXf_dyn QtQ = Q_star.transpose() * Q_star;
            MatrixXf_dyn QtX = Q_star.transpose() * X_star;
            MatrixXf_dyn XtQ = QtX.transpose();
            
            R_star = X_star.transpose() * X_star - XtQ * QtQ.llt().solve(QtX);
            
            cout << "Computed covariate-adjusted correlation matrix" << endl;
        } else {
            R_star = X_star.transpose() * X_star;
            cout << "Computed correlation matrix" << endl;
        }
        
        // Step 7: Compute effective sample size
        float n_eff = R_star.diagonal().mean();
        
        // Step 8: Normalize to correlation matrix (reuse R_star memory)
        VectorXf_dyn diag_sqrt = R_star.diagonal().array().sqrt();
        
        // Use efficient array operations instead of nested loops
        R_star = R_star.array().rowwise() / diag_sqrt.transpose().array();
        R_star = R_star.array().colwise() / diag_sqrt.array();
        
        // R_star is now the correlation matrix (R_adj)
        
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        
        cout << "LD adjustment computation completed in " << duration.count() << " ms" << endl;
        
        // Write outputs
        writeResults(output_prefix, R_star, n_eff, duration.count());
        
        return true;
    }
    
private:
    void writeResults(const string& prefix, const MatrixXf_dyn& R_adj, float n_eff, long computation_time_ms) {
        // Calculate realistic memory usage 
        size_t original_grm_memory = original_grm->rows() * original_grm->cols() * sizeof(float);
        size_t grm_subset_memory = grm_subset_indices.size() * grm_subset_indices.size() * sizeof(float);
        size_t snp_memory = R_adj.rows() * R_adj.cols() * sizeof(float);
        size_t raw_memory = raw_data.rows() * raw_data.cols() * sizeof(float);
        size_t covariate_memory = use_covariates ? covariate_data.rows() * covariate_data.cols() * sizeof(float) : 0;
        
        size_t peak_memory_mb = (original_grm_memory + grm_subset_memory + snp_memory + raw_memory + covariate_memory) / (1024 * 1024);
        
        // Write summary file
        ofstream summary_file(prefix + ".summary");
        summary_file << "Effective sample size: " << n_eff << endl;
        summary_file << "Number of matched individuals: " << matched_individuals.size() << endl;
        summary_file << "Number of SNPs: " << snp_ids.size() << endl;
        summary_file << "Heritability: " << h2 << endl;
        if (use_covariates) {
            summary_file << "Number of covariates: " << covariate_data.cols() << endl;
        }
        summary_file << "Computation time (ms): " << computation_time_ms << endl;
        summary_file << "Peak memory usage (MB): " << peak_memory_mb << endl;
        summary_file.close();
        
        // Write LD matrix
        ofstream ld_file(prefix + ".ld");
        ld_file << fixed << setprecision(6);
        for (int i = 0; i < R_adj.rows(); ++i) {
            for (int j = 0; j < R_adj.cols(); ++j) {
                if (j > 0) ld_file << " ";
                ld_file << R_adj(i, j);
            }
            ld_file << endl;
        }
        ld_file.close();
        
        // Write SNP IDs
        ofstream snp_file(prefix + ".snpids");
        for (const string& snp_id : snp_ids) {
            snp_file << snp_id << endl;
        }
        snp_file.close();
        
        cout << "Results written to:" << endl;
        cout << "  " << prefix << ".summary" << endl;
        cout << "  " << prefix << ".ld" << endl;
        cout << "  " << prefix << ".snpids" << endl;
    }
};

int main(int argc, char* argv[]) {
    string raw_file, grm_prefix, output_prefix;
    string error_weight_file, error_weight_name;
    string covariate_file, covariate_names_str;
    float h2 = 0.5f;
    int num_threads = 1;
    
    cout << "LD Adjuster v1.1.1 by Jicai Jiang" << endl;
    cout << "=================================" << endl;
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--raw" && i + 1 < argc) {
            raw_file = argv[++i];
        } else if (arg == "--grm" && i + 1 < argc) {
            grm_prefix = argv[++i];
        } else if (arg == "--h2" && i + 1 < argc) {
            h2 = stof(argv[++i]);
        } else if (arg == "--out" && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if ((arg == "--threads" || arg == "--num_threads") && i + 1 < argc) {
            num_threads = stoi(argv[++i]);
        } else if (arg == "--error_weight_file" && i + 1 < argc) {
            error_weight_file = argv[++i];
        } else if (arg == "--error_weight_name" && i + 1 < argc) {
            error_weight_name = argv[++i];
        } else if (arg == "--covariate_file" && i + 1 < argc) {
            covariate_file = argv[++i];
        } else if (arg == "--covariate_names" && i + 1 < argc) {
            covariate_names_str = argv[++i];
        } else if (arg == "--help") {
            cout << "Usage: " << argv[0] << " --raw <file> --grm <prefix> --h2 <value> --out <prefix> [options]" << endl;
            cout << "Required arguments:" << endl;
            cout << "  --raw              PLINK raw file" << endl;
            cout << "  --grm              GRM file prefix (without .grm.iid/.grm.bin)" << endl;
            cout << "  --h2               Heritability value (default: 0.5)" << endl;
            cout << "  --out              Output file prefix" << endl;
            cout << "Optional arguments:" << endl;
            cout << "  --threads N        Number of threads to use (default: 1)" << endl;
            cout << "  --num_threads N    Alias for --threads" << endl;
            cout << "  --error_weight_file CSV file with error weights (requires --error_weight_name)" << endl;
            cout << "  --error_weight_name Column name for error weights (requires --error_weight_file)" << endl;
            cout << "  --covariate_file   CSV file with covariates (requires --covariate_names)" << endl;
            cout << "  --covariate_names  Comma-separated covariate names or 'all' (requires --covariate_file)" << endl;
            cout << "  --help             Show this help message" << endl;
            return 0;
        }
    }
    
    // Validate required arguments
    if (raw_file.empty() || grm_prefix.empty() || output_prefix.empty()) {
        cerr << "Error: Missing required arguments. Use --help for usage." << endl;
        return 1;
    }
    
    // Validate paired optional arguments
    bool has_error_weight_file = !error_weight_file.empty();
    bool has_error_weight_name = !error_weight_name.empty();
    bool has_covariate_file = !covariate_file.empty();
    bool has_covariate_names = !covariate_names_str.empty();
    
    if (has_error_weight_file && !has_error_weight_name) {
        cerr << "Error: --error_weight_file requires --error_weight_name" << endl;
        return 1;
    }
    if (has_error_weight_name && !has_error_weight_file) {
        cerr << "Error: --error_weight_name requires --error_weight_file" << endl;
        return 1;
    }
    if (has_covariate_file && !has_covariate_names) {
        cerr << "Error: --covariate_file requires --covariate_names" << endl;
        return 1;
    }
    if (has_covariate_names && !has_covariate_file) {
        cerr << "Error: --covariate_names requires --covariate_file" << endl;
        return 1;
    }
    
    // Set number of threads
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
    Eigen::setNbThreads(num_threads);
#else
    Eigen::setNbThreads(num_threads);
#endif
    
    cout << "Raw file: " << raw_file << endl;
    cout << "GRM prefix: " << grm_prefix << endl;
    cout << "Heritability: " << h2 << endl;
    cout << "Output prefix: " << output_prefix << endl;
    cout << "Threads: " << num_threads << endl;
    
    if (has_error_weight_file) {
        cout << "Error weight file: " << error_weight_file << endl;
        cout << "Error weight column: " << error_weight_name << endl;
    }
    
    if (has_covariate_file) {
        cout << "Covariate file: " << covariate_file << endl;
        cout << "Covariate names: " << covariate_names_str << endl;
    }
    
    cout << endl;
    
    // Load data
    PLINKRawReader raw_reader;
    if (!raw_reader.read(raw_file)) {
        return 1;
    }
    
    GRMReader grm_reader;
    if (!grm_reader.read(grm_prefix)) {
        return 1;
    }
    
    // Load optional CSV files
    CSVReader* weight_reader = nullptr;
    CSVReader* covariate_reader = nullptr;
    vector<string> covariate_col_names;
    
    if (has_error_weight_file) {
        weight_reader = new CSVReader();
        if (!weight_reader->read(error_weight_file)) {
            delete weight_reader;
            return 1;
        }
        
        // Check if the specified column exists
        if (!weight_reader->hasColumn(error_weight_name)) {
            cerr << "Error: Column '" << error_weight_name << "' not found in error weight file" << endl;
            delete weight_reader;
            return 1;
        }
    }
    
    if (has_covariate_file) {
        covariate_reader = new CSVReader();
        if (!covariate_reader->read(covariate_file)) {
            delete weight_reader;
            delete covariate_reader;
            return 1;
        }
        
        // Parse covariate column names
        if (covariate_names_str == "all") {
            covariate_col_names = covariate_reader->column_names;
        } else {
            istringstream iss(covariate_names_str);
            string token;
            while (getline(iss, token, ',')) {
                // Trim whitespace
                token.erase(0, token.find_first_not_of(" \t\r\n"));
                token.erase(token.find_last_not_of(" \t\r\n") + 1);
                
                if (!covariate_reader->hasColumn(token)) {
                    cerr << "Error: Column '" << token << "' not found in covariate file" << endl;
                    delete weight_reader;
                    delete covariate_reader;
                    return 1;
                }
                covariate_col_names.push_back(token);
            }
        }
        
        if (covariate_col_names.empty()) {
            cerr << "Error: No valid covariate columns specified" << endl;
            delete weight_reader;
            delete covariate_reader;
            return 1;
        }
    }
    
    // Perform LD adjustment
    LDAdjuster adjuster(h2);
    if (!adjuster.loadData(raw_reader, grm_reader, weight_reader, error_weight_name, 
                          covariate_reader, covariate_col_names)) {
        delete weight_reader;
        delete covariate_reader;
        return 1;
    }
    
    if (!adjuster.computeLDAdjustment(output_prefix)) {
        delete weight_reader;
        delete covariate_reader;
        return 1;
    }
    
    // Clean up
    delete weight_reader;
    delete covariate_reader;
    
    cout << "LD adjustment completed successfully!" << endl;
    return 0;
}
