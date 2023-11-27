#include "vectocirc.h"

template <typename T>
VectoCirc<T>& VectoCirc<T>::operator=(VectoCirc<T>&& other) noexcept
{
	data_ = move(other.data_);
	start_ = other.start_;
	end_ = other.end_;
	other.start_ = other.end_ = 0;
	other.data_.clear();
	return *this;
}

template <typename T>
VectoCirc<T>::VectoCirc(VectoCirc<T>&& other) noexcept
{
	data_ = move(other.data_);
	start_ = other.start_;
	end_ = other.end_;
	other.start_ = other.end_ = 0;
	other.data_.clear();
}

template <typename T>
VectoCirc<T>::VectoCirc(VectoCirc<T>& other) noexcept
{
	data_ = other.data_;
	start_ = other.start_;
	end_ = other.end_;
}

template <typename T>
VectoCirc<T>::VectoCirc(std::vector<T>& WAVE) : start_(0), data_(WAVE), end_(data_.size())
{
	if (data_.size() == data_.capacity())
		throw std::logic_error("vector is full");
	data_.resize(data_.capacity());
}

template <typename T>
VectoCirc<T>::VectoCirc(std::vector<T>&& WAVE) : start_(0), data_(move(WAVE)), end_(data_.size())
{
	if (data_.size() == data_.capacity())
		throw std::logic_error("vector is full");
	data_.resize(data_.capacity());
}

template <typename T>
T& VectoCirc<T>:: operator[](size_t index)
{
	return data_[I(index)];
}

template <typename T>
const T& VectoCirc<T>:: operator[](size_t index) const
{
	return data_[I(index)];
}

template <typename T>
int VectoCirc<T>::frontmove() const
{
	return front_move;
}

template <typename T>
int VectoCirc<T>::endmove() const
{
	return end_move;
}

template <typename T>
void VectoCirc<T>::Strelief() const
{
	end_move = front_move = 0;
}

template <typename T>
void VectoCirc<T>::extend_front(size_t n)
{
	if (!n) return;
	start_ = (start_ - n + data_.size()) % data_.size();
	front_move -= n;
	std::fill(begin(), begin() + n, 0);
}

template <typename T>
void VectoCirc<T>::extend_back(size_t n)
{
	if (!n) return;
	end_ = (end_ + n) % data_.size();
	end_move += n;
	std::fill(end() - n, end(), 0);
}

template <typename T>
void VectoCirc<T>::shrink_front(size_t n)
{
	if (!n) return;
	start_ = (start_ + n) % data_.size();
	front_move += n;
}

template <typename T>
void VectoCirc<T>::shrink_back(size_t n)
{
	if (!n) return;
	end_ = (end_ - n + data_.size()) % data_.size();
	end_move -= n;
}

template <typename T>
void VectoCirc<T>::push_front(const VectoCirc<T>& array)
{
	size_t n = array.size();
	if (!n) return;
	extend_front(n);
	for (size_t i = 0; i < n; ++i)
		data_[I(i)] = array[i];
}

template <typename T>
void VectoCirc<T>::push_front(const std::vector<T>& array)
{
	size_t n = array.size();
	if (!n) return;
	extend_front(n);
	for (size_t i = 0; i < n; ++i)
		data_[I(i)] = array[i];
}

template <typename T>
void VectoCirc<T>::push_front(const T& element)
{
	extend_front(1);
	data_[I(0)] = element;
}

template <typename T>
void VectoCirc<T>::push_back(const VectoCirc<T>& array)
{
	size_t n = array.size();
	if (!n) return;
	extend_back(n);
	for (size_t i = size() - n; i < size(); ++i)
		data_[I(i)] = array[i - size() + n];
}

template <typename T>
void VectoCirc<T>::push_back(const std::vector<T>& array)
{
	size_t n = array.size();
	if (!n) return;
	extend_back(n);
	for (size_t i = size() - n; i < size(); ++i)
		data_[I(i)] = array[i - size() + n];
}

template <typename T>
void VectoCirc<T>::push_back(const T& element)
{
	extend_back(1);
	data_[I(size() - 1)] = element;
}

template<typename T>
T VectoCirc<T>::pop_front()
{
	T temp = data_[start_];
	start_ = (start_ + 1) % data_.size();
	++front_move;
	return temp;
}

template<typename T>
T VectoCirc<T>::pop_back()
{
	T temp = data_[(end_ - 1 + data_.size()) % data_.size()];
	end_ = (end_ - 1 + data_.size()) % data_.size();
	--end_move;
	return temp;
}

template<typename T>
void VectoCirc<T>::slide_back(const std::vector<T>& array)
{
	size_t n = array.size();
	if (n > size())
		throw std::logic_error("array size is bigger than vector size");
	shrink_front(n);
	extend_back(n);
	for (size_t i = 0; i < n; ++i)
		data_[I(i + size() - n)] = array[i];
}


template<typename T>
void VectoCirc<T>::slide_front(const std::vector<T>& array)
{
	size_t n = array.size();
	if (n > size())
		throw std::logic_error("array size is bigger than vector size");
	if (!n) return;
	shrink_back(n);
	extend_front(n);
	for (size_t i = 0; i < n; ++i)
		data_[I(i)] = array[i];

}

template <typename T>
size_t VectoCirc<T>::size() const
{
	return (end_ - start_ + data_.size()) % data_.size();
}

template <typename T>
void VectoCirc<T>::clear()
{
	start_ = end_ = 0;
	data_.clear();
}

template <typename T>
TailVectoCirc<T>::TailVectoCirc(TailVectoCirc& other) :vectorcirc(other.vectorcirc), reservetail(other.reservetail) {}

template <typename T>
TailVectoCirc<T>::TailVectoCirc(TailVectoCirc&& other) noexcept :vectorcirc(std::move(other.vectorcirc)),
reservetail(other.reservetail), virtualfront(other.virtualfront), virtualflag(other.virtualflag) { other.clear(); };

template <typename T>
TailVectoCirc<T>& TailVectoCirc<T>::operator=(TailVectoCirc&& other) noexcept
{
	vectorcirc = std::move(other.vectorcirc);
	reservetail = other.reservetail;
	virtualfront = other.virtualfront;
	virtualflag = other.virtualflag;
	other.clear();
	return *this;
}

template <typename T>
TailVectoCirc<T>::TailVectoCirc(VectoCirc<T>& array, size_t vf) :vectorcirc(array), reservetail(vf) {};

template <typename T>
TailVectoCirc<T>::TailVectoCirc(VectoCirc<T>&& array, size_t vf) :vectorcirc(std::move(array)), reservetail(vf) {};

template <typename T>
void TailVectoCirc<T>::virtual_transform(size_t rt)
{
	reservetail = rt;
	virtual_transform();
}

template <typename T>
void TailVectoCirc<T>::virtual_transform()
{
	if (virtualflag)
		return;
	if (!reservetail)
		throw std::logic_error("reservetail is 0");
	if (reservetail > vectorcirc.size())
		throw std::logic_error("reservetail is bigger than size");
	if (reservetail == vectorcirc.size())
		return;
	virtualflag = true;
	if (vectorcirc.start_ < vectorcirc.end_)
	{
		for (size_t i = 0; i < reservetail; ++i)
			vectorcirc.data_[i] = vectorcirc.data_[vectorcirc.start_ + i + vectorcirc.size() - reservetail];
		virtualfront = vectorcirc.size() - reservetail;
		vectorcirc.data_.resize(reservetail + ADD_CIR);
		vectorcirc.data_.shrink_to_fit();
		vectorcirc.start_ = 0;
		vectorcirc.end_ = reservetail;
	}
	else
	{
		std::vector<T> temp(reservetail + ADD_CIR);
		for (size_t i = 0; i < reservetail; ++i)
			temp[i] = vectorcirc[i + vectorcirc.size() - reservetail];
		virtualfront = size() - reservetail;
		vectorcirc.start_ = 0;
		vectorcirc.end_ = reservetail;
		vectorcirc.data_.swap(temp);
		temp.clear();
	}
}

template<typename T>
T& TailVectoCirc<T>::operator[](size_t index)
{
	return vectorcirc[index - virtualfront];
}

template<typename T>
const T& TailVectoCirc<T>::operator[](size_t index) const
{
	return vectorcirc[index - virtualfront];
}

template<typename T>
int TailVectoCirc<T>::frontmove() const
{
	return vectorcirc.frontmove();
}

template<typename T>
int TailVectoCirc<T>::endmove() const
{
	return vectorcirc.endmove();
}

template<typename T>
void TailVectoCirc<T>::Strelief() const
{
	vectorcirc.Strelief();
}

template<typename T>
void TailVectoCirc<T>::extend_back(size_t n)
{
	vectorcirc.extend_back(n);
}

template<typename T>
void TailVectoCirc<T>::shrink_front(size_t n)
{
	vectorcirc.shrink_front(n);
}

template<typename T>
void TailVectoCirc<T>::shrink_back(size_t n)
{
	vectorcirc.shrink_back(n);
}

template<typename T>
void TailVectoCirc<T>::push_back(const VectoCirc<T>& array)
{
	vectorcirc.push_back(array);
}

template<typename T>
void TailVectoCirc<T>::push_back(const std::vector<T>& array)
{
	vectorcirc.push_back(array);
}

template<typename T>
void TailVectoCirc<T>::push_back(const T element)
{
	vectorcirc.push_back(element);
}

template<typename T>
void TailVectoCirc<T>::slide_back(const std::vector<T>& array)
{
	vectorcirc.slide_back(array);
}

template<typename T>
T TailVectoCirc<T>::pop_front()
{
	return vectorcirc.pop_front();
}

template<typename T>
T TailVectoCirc<T>::pop_back()
{
	return vectorcirc.pop_back();
}

template<typename T>
size_t TailVectoCirc<T>::size() const
{
	return vectorcirc.size() + virtualfront;
}

template<typename T>
void TailVectoCirc<T>::clear()
{
	vectorcirc.clear();
}