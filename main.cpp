#include "common_functions.hpp"

using namespace std;

// @return	the index the column that holds times
size_t findDAColumn(const vector<string> & column_names) {
    const regex DA_column("DA:.*");
    size_t column_no = 0;
    for (const string & name : column_names) {
        if (regex_match(name, DA_column)) {
            return column_no;
        }
        column_no++;
    }

    throw runtime_error("DA_column not found.\n");
}

// @return	a name of a component from a column name
string obtainName(const string & column_name) {
    if (regex_match(column_name, regex("DV:.*")))
        return column_name.substr(3, column_name.size() - 3);
    if (regex_match(column_name, regex("TR:.*:Inhibitors")))
        return column_name.substr(3, column_name.size() - 3 - strlen("Inhibitors") - 1);
    if (regex_match(column_name, regex("TR:.*:Stimuli")))
        return column_name.substr(3, column_name.size() - 3 - strlen("Stimuli") - 1);
    if (regex_match(column_name, regex("TR:.*(i|a)")))
        return column_name.substr(3, column_name.size() - 4);
    return "";
}

// @return	true iff the column holds a component
bool isComponent(const string & column_name) {
    return !obtainName(column_name).empty();
}

// @return	true if the component is included in some measurement
bool isMeasured(const string & column_name) {
    const regex component_expr("DV:.*");
    return regex_match(column_name, component_expr);
}

// A structure that holds results for a component
struct ComponentData {
    size_t column_no; ///< Index of the column holding its values.
    string name; ///< Name of the component.
    bool measured; ///< True iff the component is not either stimulated or inhibited.
    set<string> inhibitors; ///< Names of the species that inhibit this node.
    set<string> activators; ///< Names of the species that activate this node.
};

// @return	a vector of components holding the data relevant to the MIDAS file
vector<ComponentData> getComponenets(const vector<string> & column_names) {
    vector<ComponentData>  components;
    size_t column_no = 0;
    for (const string & column : column_names) {
        if (isComponent(column)) {
            components.push_back({column_no, obtainName(column), isMeasured(column), set<string>(), set<string>()});
        }
        column_no++;
    }
    return components;
}

// @return	a 2D vector (first columns then row) containing the data from the MIDAS file
vector<vector<string> > getData(fstream & input_file) {
    vector<vector<string > > data;

    string line;
    while(getline(input_file, line)) {
        vector<string> line_data;
        boost::split(line_data, line, boost::is_any_of(","));
        data.push_back(line_data);
    }

    return data;
}

// @return	the time value (expecte to be integer) in the given row
inline size_t getTime(const vector<string> & row, const size_t DA_column) {
	try {
		return boost::lexical_cast<size_t>(row[DA_column]);
	}
	catch (...) {
		throw invalid_argument("Non-integral time in the input file.");
	}
}

// @return	a set of all the timepoints
set<size_t> getTimes(const vector<vector<string> > & data, const size_t DA_column) {
    set<size_t> times;

    for (const vector<string> line : data) {
        times.insert(getTime(line, DA_column));
    }

    return times;
}

// @return	true iff the given row occurs in the requested timepoint
inline bool isInTime(const vector<string> row, const size_t DA_column, const size_t time) {
	return (getTime(row, DA_column) == time);
}

// @return	true iff the source and targe match in the experiment setup (inhibitions and activations)
bool matchesSetting(const vector<ComponentData> & components, const vector<string> & source, const vector<string> target) {
    for (const size_t comp_no : cscope(components)) {
        if (components[comp_no].measured)
            continue;

        if (source[comp_no].compare(source[comp_no]) != 0)
            return false;
    }
    return true;
}

// @return	pairs of rows that represent consecutive measurement
vector<pair<vector<string>, vector<string> > > 
getTimesteps(const vector<ComponentData> & components, const size_t DA_column, const vector<vector<string> > & data) {
    vector<pair<vector<string>, vector<string> > > timesteps; // Measurements are grouped by a timepoint, not by the experiment!

    const set<size_t> times = getTimes(data, DA_column);

    for (const vector<string> & source_line : data) {
        const size_t source_time = getTime(source_line, DA_column);
        auto time_it = times.find(source_time);
        if(++time_it == times.end())
            continue;

		// Get the successor for the current line and push the two in the vector if there is such
        for (const vector<string> & target_line : data) {
            if (!isInTime(target_line, DA_column, *time_it))
                continue;

            if (matchesSetting(components, source_line, target_line)) {
                timesteps.push_back(make_pair(source_line, target_line));
                break;
            }
        }
    }

    return timesteps;
}

// @return	vector of timesteps that no longer hold timing information
vector<pair<vector<double>, vector<double> > > 
convert(const vector<ComponentData> & components, const vector<pair<vector<string>, vector<string> > > & timesteps) {
    vector<pair<vector<double>, vector<double> > > converted;

    for (const pair<vector<string>, vector<string> >  & timestep : timesteps) {
        vector<double> first;
        vector<double> second;
        for (const ComponentData & component : components) {
            first.push_back(boost::lexical_cast<double>(timestep.first[component.column_no]));
            second.push_back(boost::lexical_cast<double>(timestep.second[component.column_no]));
        }
        converted.push_back(make_pair(first, second));
    }

    return converted;
}

// @return	vector of timesteps that has been normalised
// normalisation means that all the ranges were narrowed to [0,1] and the values has been adjusted accordingly
vector<pair<vector<double>, vector<double> > > normalize(const vector<pair<vector<double>, vector<double> > > & timesteps) {
    vector<pair<vector<double>, vector<double> > > normalized;

    if (timesteps.empty())
        throw runtime_error("Nothing to normalize.");
    const size_t comp_count = timesteps.front().first.size();

	// Compute boundaries
    vector<double> minimals(comp_count, numeric_limits<double>::max());
    vector<double> maximals(comp_count, 0.0);
    for (const pair<vector<double>, vector<double> >  & timestep : timesteps) {
        for (const size_t comp_no : crange(comp_count)) {
            minimals[comp_no] = min(minimals[comp_no], min(timestep.first[comp_no], timestep.second[comp_no]));
            maximals[comp_no] = max(maximals[comp_no], max(timestep.first[comp_no], timestep.second[comp_no]));
        }
    }

	// Compute differences between boundaries
    vector<double> difference;
    for (const size_t comp_no : crange(comp_count)) {
        difference.push_back(maximals[comp_no] - minimals[comp_no]);
    }

    for (const pair<vector<double>, vector<double> >  & timestep : timesteps) {
        vector<double> first, second;
        for (const size_t comp_no : crange(comp_count)) {
            first.push_back((timestep.first[comp_no] - minimals[comp_no]) / difference[comp_no]);
            second.push_back((timestep.second[comp_no] - minimals[comp_no]) / difference[comp_no]);
        }
        normalized.push_back(make_pair(first, second));
    }

    return normalized;
}

struct Interval {
    const vector<double> first_;
    const vector<double> second_;
    vector<double> changes_;
    vector<double> productions_;

    NO_COPY_SHORT(Interval)
    DEFAULT_MOVE(Interval)

    Interval(const vector<double> & first, const vector<double> & second) : first_(first), second_(second) {
        double avg_production = 0;
        size_t production_cases = 0;

        // Compute change
        for (const size_t val_no : cscope(first)) {
            const double change = second[val_no] - first[val_no];
            changes_.push_back(change);
            if ((change) > 0.) {
                avg_production += change;
                production_cases++;
            }
        }
        avg_production /= production_cases;

        // Compute production factor - positive if bigger than average, negative otherwise
        for (const double change : changes_) {
            if (change <= 0.)
                productions_.push_back(0.);
            else
                productions_.push_back(change - avg_production);
        }
    }
};

// Compute time intervals between measurements
vector<Interval> makeIntervals(const vector<pair<vector<double>, vector<double> > > & timesteps) {
    vector<Interval> intervals;

    for (const pair<vector<double>, vector<double> >  & timestep : timesteps) {
        intervals.push_back(Interval(timestep.first, timestep.second));
    }

    return intervals;
}

// Standard logistic curve (expects x to be in [0,1])
double inline logistic(const double x) {
    return 1 / (1 + exp(-(x*12-6)));
}

// Average of x and y weighted by the diversity between the two
double inline logistic(const double x, const double y) {
    return logistic((max(x,y) + x*y)/2);
}

// Average of x, y,z  weighted by the diversity between the three
double inline logistic(const double x, const double y, const double z) {
    return logistic((max(x,max(z,y)) + x*y + x*z + y*z - y*x*z)/3);
}

// Compute regulators (positive or negative) of all the components
vector<ComponentData> computeRegulators(const vector<ComponentData> & components, vector<Interval> & intervals, bool inhibition) {
    vector<ComponentData> regulated = components;

	// Counting values for output
    const size_t regulated_count = count_if(components.begin(), components.end(), [](const ComponentData & comp){return comp.measured;}) ;
    // The second part is choose three from components.size() with repetitions
    const size_t total = regulated_count * ((components.size() * (components.size() + 1) * (components.size() + 2)) / 3);
    size_t current = inhibition ? 0 : total / 2;

	// Compute for all the components
    for (const size_t target_no : cscope(components)) {
        map<set<string>,  double> mutual_inf;
        if (!components[target_no].measured)
            continue;

        // Test all combinations up to three regulators (combination of multiple selfs == just single one).
        for (const size_t source_1_no : cscope(components)) {
            for (const size_t source_2_no : crange(source_1_no, components.size())) {
                for (const size_t source_3_no : crange(source_2_no, components.size())) {
                    cout << "Testing: " << current++ << "/" << total << "     \r";
                    double effect = 0.;
                    for (const Interval & interval : intervals) {
                        double factor = inhibition ? -1 * interval.changes_[target_no] : interval.productions_[target_no];
                        effect += factor * logistic(interval.first_[source_1_no],interval.first_[source_2_no],interval.first_[source_3_no]);
                    }
                    mutual_inf.insert(
						make_pair(set<string>({components[source_1_no].name, components[source_2_no].name, components[source_3_no].name}), effect));
                }
            }
        }

        auto max_it = max_element(mutual_inf.begin(), mutual_inf.end(), [](const pair<set<string>, double> & A, const pair<set<string>, double> & B){
            return A.second < B.second;
        });

        // Pick the set of regulstors with the biggest positive influence - if there is none, just take empty.
        set<string> & sources = inhibition ? regulated[target_no].inhibitors : regulated[target_no].activators;
        const double PENALTY = -0.00;
        sources = (max_it->second  + PENALTY * max_it->first.size()) > 3*PENALTY ? max_it->first : set<string>();
    }

    return regulated;
}


void output(const vector<ComponentData> & components, const string & file_name) {
	string new_filename = file_name.substr(0, file_name.length() - 3) + "sif";
	ofstream fout(new_filename, ios::out);

    for (const ComponentData & component: components) {
        for (const string & source : component.activators) {
            fout << source << " 1 " << component.name << endl;
        }
        for (const string & source : component.inhibitors) {
            fout << source << " -1 " << component.name << endl;
        }
    }
}

// The main function expects a csv file in the MIDAS format
// The function outpus a .sif graph file in the same folder and with the same name
int main(int argc, char* argv[]) {
	// Get input file
    string filename(argv[1]);
    fstream input_file(filename, ios::in);
    if (!input_file)
        throw invalid_argument("Wrong filename \"" + filename + "\".\n");

	// Read column names
    string names_line;
    getline(input_file, names_line);
    vector<string> column_names;
    boost::split(column_names, names_line, boost::is_any_of(","));

	// Obtain data and normalize them
    size_t DA_colum = findDAColumn(column_names);
    vector<ComponentData> components = getComponenets(column_names);
    vector<vector<string> > data = getData(input_file);
    vector<pair<vector<string>, vector<string> > > timepoints = getTimesteps(components, DA_colum, data);
    vector<pair<vector<double>, vector<double> > > converted = convert(components, timepoints);
    vector<pair<vector<double>, vector<double> > > normalized = normalize(converted);
    vector<Interval> intervals = makeIntervals(normalized);

	// Add positive and negative regulators
    components = computeRegulators(components, intervals, true);
    components = computeRegulators(components, intervals, false);
	
	// Output
	output(components, filename);

    return 0;
}

