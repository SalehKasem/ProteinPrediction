#pragma once

class Node {
protected:
#pragma region Protected vars
    /*
     * protein index
     */
    long m_protein;
    /*
     * fragment index
     */
    int m_index;
    /*
     * hamming distance
     */
    int hamming;
    /*
     * contex hamming  distance
     */
    int contex;
    /*
     * amino similarity with neighbor node (as neighbor)
     */
    double aminoSeq;
    /*
     * rmsd with neighbor node (as neighbor)
     */
    double rmsd;
    /*
     * mean rmsd with neighbor node (as neighbor) -Linear Regression
     */
    double meanRmsdLinear, meanRmsdLinear2;
    /*
     * mean rmsd with neighbor node (as neighbor) - Logistic Regression
     */
    double meanRmsdLogistic, meanRmsdLogistic2;
    /*
     * mean rmsd with neighbor node (as neighbor) - Monotonic Regression
     */
    double meanRmsdMonotonic;
    /*
     * weight with neighbor node calculated by weight function
     */
    double weight;
    /*
     * type of similarity sequence
     */
    int type;
#pragma endregion Protected vars
public:
    Node();

    Node(long proteind, int index);

    Node(long proteind, int index, int hamming, int contex, double rmsd, double meanRmsdLinear, double meanRmsdLogistic,
         int type, double aminoSqe);

    ~Node();

#pragma region Getter Setter

    /* getter and setter section Start*/
    long getProteinIndex();

    void setProteinIndex(long m_currentProtein);

    int getFragmentIndex();

    void setFragmentIndex(int m_index);

    double getMeanRmsdLinear();

    void setMeanRmsdLinear(double m_meanRmsd);

    double getMeanRmsdLogistic();

    void setMeanRmsdLogistic(double m_meanRmsd);

    void setWeight(double weight);

    double getWeight();

    int getHamming();

    void setHamming(int hamming);

    int getContex();

    void setContex(int contex);

    double getRmsd();

    void setRmsd(double rmsd);

    double getAminoSeq();

    void setAminoSeq(double aminoSeq);

    int getType();

    void setType(int type);

    double getMeanRmsdMonotonic();

    void setMeanRmsdMonotonic(double meanRmsdMonotonic);

    double getMeanRmsdLinear2();

    void setMeanRmsdLinear2(double meanRmsdLinear2);

    double getMeanRmsdLogistic2();

    void setMeanRmsdLogistic2(double meanRmsdLogistic2);
    /* getter and setter section END*/
#pragma endregion Getter Setter
};

