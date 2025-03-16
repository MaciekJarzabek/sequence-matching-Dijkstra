#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <random>
#include <algorithm>
#include <queue>
#include <limits>

#define inf 1000
using namespace std;

struct Node {
    int vertex;
    int cost;

    bool operator>(const Node& other) const {
        return cost > other.cost;
    }
};

int len_olig;

string generateRandomSequence(int length) {
    string sequence = "";
    string nucleotides = "ACGT";
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> distribution(0, nucleotides.size() - 1);

    for (int i = 0; i < length; ++i) {
        sequence += nucleotides[distribution(generator)];
    }

    return sequence;
}

vector<string> splitSequence(const string& sequence, int subLength) {
    vector<string> subSequences;
    for (int i = 0; i <= sequence.length() - subLength; ++i) {
        subSequences.push_back(sequence.substr(i, subLength));
    }
    return subSequences;
}

vector<string> removeDuplicates(vector<string>& subSequences) {
    set<string> uniqueSet(subSequences.begin(), subSequences.end());
    subSequences.assign(uniqueSet.begin(), uniqueSet.end());
    return subSequences;
}

int calculateMisplaced(const string& seq1, const string& seq2) {
    int misplaced = 0;
    for (int i = 0; i < seq1.length(); ++i) {
        if (seq1[i] != seq2[i]) {
            ++misplaced;
        }
    }
    return misplaced;
}

map<string, vector<pair<string, int>>> createHybridizationGraph(const vector<string>& ver_info, int len_olig) {
    map<string, vector<pair<string, int>>> graph;
    int vertex_amount = ver_info.size();
    vector<vector<int>> g_matrix;

    g_matrix.resize(vertex_amount, vector<int>(vertex_amount, inf));
    for (int i = 0; i < vertex_amount; i++) {
        for (int j = 0; j < vertex_amount; j++) {
            if (i != j) {
                string first = ver_info[i].substr(1);
                string second = ver_info[j].substr(0, len_olig - 1);

                if (first == second) {
                    g_matrix[i][j] = 1;
                }

                first = ver_info[i].substr(2);
                second = ver_info[j].substr(0, len_olig - 2);

                if (first == second) {
                    g_matrix[i][j] = 2;
                }
            }
        }
    }

    for (int i = 0; i < vertex_amount; ++i) {
        for (int j = 0; j < vertex_amount; ++j) {
            if (i != j && g_matrix[i][j] > 0) {
                string node1 = ver_info[i];
                string node2 = ver_info[j];

                int misplaced = g_matrix[i][j];
                graph[node1].push_back({ node2, misplaced });
            }
        }
    }

    return graph;
}

void printGraph(const map<string, vector<pair<string, int>>>& graph) {
    int misplaced1 = 0;
    int misplaced2 = 0;

    for (const auto& node : graph) {
        cout << "Podsekwencja: " << node.first << '\n';
        for (const auto& edge : node.second) {
            if (edge.second == 1 || edge.second == 2) {
                cout << "  Laczy sie z: " << edge.first << ", misplaced: " << edge.second << '\n';
                if (edge.second == 1) {
                    ++misplaced1;
                }
                else if (edge.second == 2) {
                    ++misplaced2;
                }
            }
        }
    }

    cout << "Liczba misplaced 1: " << misplaced1 << '\n';
    cout << "Liczba misplaced 2: " << misplaced2 << '\n';
}

string reconstructSequence(const map<string, vector<pair<string, int>>>& graph, const vector<string>& subSequences, int originalLength, int positiveErrors) {
    string reconstructedSequence;
    set<string> visited;

    for (const string& startNode : subSequences) {
        if (visited.find(startNode) != visited.end()) {
            continue;
        }

        visited.insert(startNode);
        reconstructedSequence += startNode;

        while (true) {
            const string& currentSequence = reconstructedSequence.substr(reconstructedSequence.length() - len_olig + 1);

            auto it = graph.find(currentSequence);
            if (it != graph.end() && !it->second.empty()) {
                const vector<pair<string, int>>& edges = it->second;

                vector<pair<string, int>> sortedEdges(edges.begin(), edges.end());
                sort(sortedEdges.begin(), sortedEdges.end(), [](const auto& a, const auto& b) {
                    return a.second < b.second;
                    });

                bool found = false;
                for (const auto& edge : sortedEdges) {
                    const string& nextNode = edge.first;
                    if (visited.find(nextNode) == visited.end()) {
                        int overlap = edge.second;
                        reconstructedSequence += nextNode.substr(overlap);
                        visited.insert(nextNode);
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    break;
                }
            }
            else {
                break;
            }
        }
    }


    reconstructedSequence.resize(originalLength, ' ');

    return reconstructedSequence;
}


int calculateLevenshteinDistance(const string& seq1, const string& seq2) {
    int m = seq1.length();
    int n = seq2.length();

    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (i == 0) {
                dp[i][j] = j;
            }
            else if (j == 0) {
                dp[i][j] = i;
            }
            else if (seq1[i - 1] == seq2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            }
            else {
                dp[i][j] = 1 + min({ dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1] });
            }
        }
    }

    return dp[m][n];
}
int calculateDifferences(const string& seq1, const string& seq2) {
    int differences = 0;
    int length = min(seq1.length(), seq2.length());

    for (int i = 0; i < length; ++i) {
        if (seq1[i] != seq2[i]) {
            ++differences;
        }
    }

    return differences;
}
double calculateCoverage(const string& originalSequence, const string& reconstructedSequence) {
    int differences = calculateDifferences(originalSequence, reconstructedSequence);
    int totalLength = originalSequence.length();

    double coveragePercentage = (1.0 - static_cast<double>(differences) / totalLength) * 100.0;
    return coveragePercentage;
}

double compareSequences(const string& originalSequence, const string& reconstructedSequence) {
    int levenshteinDistance = calculateLevenshteinDistance(originalSequence, reconstructedSequence);
    int maxLength = max(originalSequence.length(), reconstructedSequence.length());

    double matchingPercentage = 100.0 - (static_cast<double>(levenshteinDistance) / maxLength * 100.0);
    return matchingPercentage;
}

std::vector<int> Dijkstras_algorithm(const vector<vector<int>>& matrix, int start, int end) {
    int m = matrix.size();
    vector<int> d(m, numeric_limits<int>::max());
    vector<int> p(m, -1);
    vector<bool> visited(m, false);

    d[start] = 0;

    for (int count = 0; count < m - 1; ++count) {
        int u = -1;
        for (int i = 0; i < m; ++i) {
            if (!visited[i] && (u == -1 || d[i] < d[u]))
                u = i;
        }

        visited[u] = true;

        for (int v = 0; v < m; ++v) {
            if (!visited[v] && matrix[u][v] != 0 && d[v] > d[u] + matrix[u][v]) {
                d[v] = d[u] + matrix[u][v];
                p[v] = u;
            }
        }
    }

    vector<int> result;
    for (int x = end; x != -1; x = p[x]) {
        result.push_back(x);
    }
    reverse(result.begin(), result.end());

    return result;
}

int main() {
    int subLength;
    cout << "Wprowadz dlugosc podsekwencji: ";
    cin >> subLength;
    //subLength = 12;
    len_olig = subLength;
    int originalLength = 650;
    string sequence = generateRandomSequence(originalLength);
    //string sequence = "";
    cout << "Wygenerowana sekwencja: " << sequence << endl;

    vector<string> subSequences = splitSequence(sequence, subLength);

    int positiveErrors;
    int negativeErrors;
    //negativeErrors = (originalLength*2)/100;
    //positiveErrors = (originalLength*2)/100;
    cout << "Wprowadz ilosc bledow negatywnych: ";
    cin >> negativeErrors;


    negativeErrors = min(negativeErrors, static_cast<int>(subSequences.size()));

    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> distribution(0, subSequences.size() - 1);

    for (int i = 0; i < negativeErrors; ++i) {
        if (!subSequences.empty()) {
            int index = distribution(generator);
            cout << "Usunieto podsekwencje: " << subSequences[index] << endl;
            subSequences.erase(subSequences.begin() + index);
        }
    }

    cout << "Wprowadz ilosc bledow pozytywnych: ";
    cin >> positiveErrors;

    for (int i = 0; i < positiveErrors; ++i) {
        string newSubSequence = generateRandomSequence(subLength);
        subSequences.push_back(newSubSequence);
        cout << "Dodano podsekwencje: " << newSubSequence << endl;
    }
    subSequences = removeDuplicates(subSequences);

    vector<vector<int>> g_matrix(subSequences.size(), vector<int>(subSequences.size(), 0));
    for (int i = 0; i < subSequences.size(); ++i) {
        for (int j = 0; j < subSequences.size(); ++j) {
            if (i != j) {
                int misplaced = calculateMisplaced(subSequences[i], subSequences[j]);
                g_matrix[i][j] = misplaced;
            }
        }
    }

    map<string, vector<pair<string, int>>> graph = createHybridizationGraph(subSequences, len_olig);

    printGraph(graph);
    subSequences = removeDuplicates(subSequences);
    int startVertex = 0;
    int endVertex = subSequences.size() - 1;

    vector<int> dijkstraPath = Dijkstras_algorithm(g_matrix, startVertex, endVertex);

    cout << "\nDijkstra's Shortest Path: ";
    for (int vertex : dijkstraPath) {
        cout << vertex << " ";
    }

    string reconstructedSequence = reconstructSequence(graph, subSequences, originalLength, positiveErrors);
    cout << "\nReconstructed Sequence: " << reconstructedSequence << endl;

    double matchingPercentage = compareSequences(sequence, reconstructedSequence);
    cout << "Matching Percentage (Levenshtein): " << matchingPercentage << "%" << endl;

    double coveragePercentage = calculateCoverage(sequence, reconstructedSequence);
    cout << "Coverage Percentage: " << coveragePercentage << "%" << endl;

    return 0;
}