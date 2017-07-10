#ifndef REOSOLVER_SPACEMESH_H
#define REOSOLVER_SPACEMESH_H

#include <vector>

struct Cell
{
	double density, pressure, velocity;
};

class SpaceMesh: public std::vector<Cell>
{

public:

	SpaceMesh(unsigned long n = 0): std::vector<Cell>(n) // or n+4 if making virtual cells.
	// making virtual cells is not required cuz of reloading [] operator
	{
		for(unsigned long i = 0; i < n; ++i)
		{
			(*this)[i].density = 0.0;
			(*this)[i].velocity = 0.0;
			(*this)[i].pressure = 0.0;
		}
	}

	Cell& operator[](int i)
	{
		//return (*this).std::vector<Cell>::operator[](i + 2);
		if(i == -1)
			return (*this).std::vector<Cell>::operator[](0);
		if(i == -2)
			return (*this).std::vector<Cell>::operator[](1);
		if(i == -3)
			return (*this).std::vector<Cell>::operator[](2);
		if(i == -4)
			return (*this).std::vector<Cell>::operator[](3);
		if(i == this->size())
			return (*this).std::vector<Cell>::operator[](static_cast<int>(this->size() - 1));
		if(i == this->size() + 1)
			return (*this).std::vector<Cell>::operator[](static_cast<int>(this->size() - 2));
		if(i == this->size() + 2)
			return (*this).std::vector<Cell>::operator[](static_cast<int>(this->size() - 3));
		if(i == this->size() + 3)
			return (*this).std::vector<Cell>::operator[](static_cast<int>(this->size() - 4));
		// this is equal to copying cells
		return (*this).std::vector<Cell>::operator[](i);
	}

	Cell const& operator[](int i) const
	{
		//return (*this).std::vector<Cell>::operator[](i + 2);
		if(i == -1)
			return (*this).std::vector<Cell>::operator[](0);
		if(i == -2)
			return (*this).std::vector<Cell>::operator[](1);
		if(i == -3)
			return (*this).std::vector<Cell>::operator[](2);
		if(i == this->size())
			return (*this).std::vector<Cell>::operator[](static_cast<int>(this->size() - 1));
		if(i == this->size() + 1)
			return (*this).std::vector<Cell>::operator[](static_cast<int>(this->size() - 2));
		if(i == this->size() + 2)
			return (*this).std::vector<Cell>::operator[](static_cast<int>(this->size() - 3));
		return (*this).std::vector<Cell>::operator[](i);
	}

/*		void CopyCell(int source, int destination)
		{
			(*this)[destination] = (*this)[source];
		} */


	friend const SpaceMesh operator-(SpaceMesh const& left, SpaceMesh const& right);


};

const SpaceMesh operator-(SpaceMesh const& left, SpaceMesh const& right);




#endif //REOSOLVER_SPACEMESH_H
