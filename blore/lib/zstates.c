#include <iostream>
#include <string>
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <vector>
using namespace std;


void print2Darray (int arr[], int rows, int columns) {
    for (int i=0; i<rows; ++i){
        for (int j=0; j<columns; ++j) {
            cout << arr[i*columns + j] << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
}

void print2Dvector (vector<int>* vec, int columns) {
    int rows = vec->size() / columns;
    cout << "Size of vectors: " << vec->size() << '\n';
    for (int i=0; i<rows; ++i){
        for (int j=0; j<columns; ++j) {
            cout << (*vec)[i*columns + j] << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
}

bool find_in_array (int x, int *a, int n) {
    for (int i=0; i<n; i++) {
        if (a[i] == x) {
            return true;
        }
    }
    return false;
}

void append_new_states (int* newz, int* lead, vector<int>* absent, int norm, int &nstate) {
    int newnorm = norm + 1;
    int new_nstate = nstate + absent->size();
    cout << "Adding states after: " << nstate << '\n';
    cout << "Absent values are:" << '\n';
    print2Dvector(absent, absent->size());

    int abscount = 0;
    for (int i=nstate; i<new_nstate; i++) {
        for (int j=0; j<norm; j++) {
            newz[i*newnorm + j] = lead[j];
        }
        cout << "Appending new z at: " << (i+1)*newnorm - 1 << " with value " << (*absent)[abscount] << '\n';
        newz[(i+1)*newnorm - 1] = (*absent)[abscount];
        abscount++;
    }
    nstate = new_nstate;
}


void create_zstates(int nstate, int norm, int nsnps, int *oldz, int *newz) {

    int newnorm = norm + 1;
    std::vector<int> absent;
    std::vector<int> present;

    bool exists;
    int tmpval;
    int tmphold;
    int statecount;
    int samecount;

    int *zlead = new int[norm];

    // Initialize with the the first lead.
    statecount = 0;
    for (int i=0; i<norm; i++) {
        zlead[i] = oldz[i];
    }

    // Get the absent snps
    for (int k=0; k<nsnps; k++) {
        exists = find_in_array(k, zlead, norm);
        if (!exists) {
            absent.push_back(k);
        }
    }

    // Get the arrays from this lead, and put it to newz
    append_new_states(newz, zlead, &absent, norm, statecount);

    // for the following leads
    for (int i=1; i<nstate; i++) {

        present.clear();
        absent.clear();

        // create the mlead
        for (int j=0; j<norm; j++) {
            zlead[j] = oldz[i*norm + j];
            present.push_back(zlead[j]);
        }

        // iterate over previous mleads to find those with (k-1) common elements
        for (int j=0; j<i; j++) {
            samecount = 0;
            for (int k=0; k<norm; k++) {
                tmpval = oldz[j*norm + k];
                exists = find_in_array(tmpval, zlead, norm);
                if (exists) {
                    samecount ++;
                } else {
                    tmphold = tmpval;
                }
            }
            if (samecount == norm - 1) { // then we have 2 leads with (k-1) similar elements
                present.push_back(tmphold);
            }
        }

        // create the absent vector
        for (int j=0; j<nsnps; j++) {
            exists = std::find(present.begin(), present.end(), j) != present.end();
            if (!exists) {
                absent.push_back(j);
            }
        }

        // create the new zstates
        append_new_states(newz, zlead, &absent, norm, statecount);
    }
    nstate = statecount;
    norm   = newnorm;
    delete [] zlead;
}

int main() {
    int nlead = 3;
    int norm = 3;
    int nsnps = 10;
    int oldz[] = {1,2,3,1,2,4,3,6,8};
    int *newz = new int[nlead * (nsnps - norm) * (norm + 1)];
    create_zstates(nlead, norm, nsnps, oldz, newz);
    print2Darray(newz, (nlead * (nsnps - norm)), norm+1);
}

// http://stackoverflow.com/questions/127704/algorithm-to-return-all-combinations-of-k-elements-from-n
