#include "NodePCN.h"


using namespace std;

vector<Node> NodePCN::getNeighbors() {
	return neighbors;
}
void NodePCN::setNeighbors(vector<Node> neighbors) {
	this->neighbors = neighbors;
}

NodePCN::NodePCN()
{
}

/**
 * constructor
 * @param protein
 * @param index
 * call super class (node) costructor
 */
NodePCN::NodePCN(long protein, int index) :Node(protein, index) {

}
/**
 * copy constructor
 * @param node
 * call constructor with protein index and , fragment index then cast node to node pcn and get neighbors 
 */
NodePCN::NodePCN(Node node) {
	NodePCN(node.getProteinIndex(), node.getFragmentIndex());
	neighbors = ((NodePCN)node).getNeighbors();
}

NodePCN::~NodePCN()
{
}