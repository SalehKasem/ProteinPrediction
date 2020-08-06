#include "Node.h"


long Node::getProteinIndex() {
    return this->m_protein;
}

void Node::setProteinIndex(long m_currentProtein) {
    this->m_protein = m_currentProtein;
}

int Node::getFragmentIndex() {
    return this->m_index;
}

void Node::setFragmentIndex(int m_index) {
    this->m_index = m_index;
}

double Node::getMeanRmsdLinear() {
    return this->meanRmsdLinear;
}

void Node::setMeanRmsdLinear(double m_meanRmsd) {
    this->meanRmsdLinear = m_meanRmsd;
}

double Node::getMeanRmsdLogistic() {
    return this->meanRmsdLogistic;
}

void Node::setMeanRmsdLogistic(double m_meanRmsd) {
    this->meanRmsdLogistic = m_meanRmsd;
}

void Node::setWeight(double weight) {
    this->weight = weight;
}

double Node::getWeight() {
    return this->weight;
}

int Node::getHamming() {
    return this->hamming;
}

void Node::setHamming(int hamming) {
    this->hamming = hamming;
}

int Node::getContex() {
    return this->contex;
}

void Node::setContex(int contex) {
    this->contex = contex;
}

double Node::getRmsd() {
    return this->rmsd;
}

void Node::setRmsd(double rmsd) {
    this->rmsd = rmsd;
}

double Node::getAminoSeq() {
    return this->aminoSeq;
}

void Node::setAminoSeq(double aminoSeq) {
    this->aminoSeq = aminoSeq;
}

int Node::getType() {
    return this->type;
}

void Node::setType(int type) {
    this->type = type;
}

double Node::getMeanRmsdMonotonic() {
    return this->meanRmsdMonotonic;
}

void Node::setMeanRmsdMonotonic(double meanRmsdMonotonic) {
    this->meanRmsdMonotonic = meanRmsdMonotonic;
}

double Node::getMeanRmsdLinear2() {
    return this->meanRmsdLinear2;
}

void Node::setMeanRmsdLinear2(double meanRmsdLinear2) {
    this->meanRmsdLinear2 = meanRmsdLinear2;
}

double Node::getMeanRmsdLogistic2() {
    return this->meanRmsdLogistic2;
}

void Node::setMeanRmsdLogistic2(double meanRmsdLogistic2) {
    this->meanRmsdLogistic2 = meanRmsdLogistic2;
}


Node::Node() {
}

/***
	 * sceond constructor
	 * @param proteind - protein index
	 * @param index - fragment index
	 */
Node::Node(long proteind, int index) {
    this->m_protein = proteind;
    this->m_index = index;
    this->setRmsd(0);
    this->weight = 0;
}

/***
 * sceond constructor
 * @param proteind - protein index
 * @param index - fragment index
 */
Node::Node(long proteind, int index, int hamming, int contex, double rmsd, double meanRmsdLinear,
           double meanRmsdLogistic, int type, double aminoSqe) {
    this->m_protein = proteind;
    this->m_index = index;
    this->setRmsd(rmsd);
    this->meanRmsdLinear = meanRmsdLinear;
    this->meanRmsdLogistic = meanRmsdLogistic;
    this->setHamming(hamming);
    this->setContex(contex);
    this->setAminoSeq(aminoSqe);
    this->setType(type);
    this->weight = 0;
}

Node::~Node() {
}