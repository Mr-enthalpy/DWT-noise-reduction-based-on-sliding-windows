#pragma once

#include <vector>
#include <compare>
#include <exception>
#include <stdexcept>

constexpr int ADD_CIR = 3;

template <typename T>
class VectoCirc;

template <typename T>
class TailVectoCirc;

template <typename T>
class VectoCirc
{ 
protected:
	std::vector<T> data_;
	size_t start_, end_;
	constexpr size_t I(size_t i)
	{
		return (start_ + i) % data_.size();
	}
	mutable long long front_move = 0, end_move = 0;
	friend class TailVectoCirc<T>;
public:
	class iterator
	{
	private:
		long long index_;
		VectoCirc<T>* outclass_;
		void check_range(long long idex) const
		{
			if (idex < 0 || idex > outclass_->size())
				throw std::out_of_range("out of range");
		}
	public:
		using iterator_category = std::random_access_iterator_tag;
		using value_type = T;
		using difference_type = std::ptrdiff_t;
		using pointer = T*;
		using reference = T&;
		iterator(VectoCirc<T>& parent, size_t index) : outclass_(&parent), index_(index) {}
		iterator(const iterator& other) : outclass_(other.outclass_), index_(other.index_) {}
		iterator(iterator&& other) noexcept : outclass_(std::move(other.outclass_)), index_(other.index_)
		{
			other.index_ = 0;
			other.outclass_ = nullptr;
		}
		reference operator*()
		{
			return (*outclass_)[index_];
		}
		iterator& operator++()
		{
			++index_;
			check_range(index_);
			return *this;
		}
		iterator operator++(int)
		{
			iterator tmp = *this;
			++index_;
			check_range(index_);
			return tmp;
		}
		iterator& operator--()
		{
			--index_;
			check_range(index_);
			return *this;
		}
		iterator operator--(int)
		{
			iterator tmp = *this;
			--index_;
			check_range(index_);
			return tmp;
		}
		iterator operator+(difference_type n) const
		{
			difference_type reindex = index_ + n;
			check_range(reindex);
			return iterator(*outclass_, reindex);
		}
		iterator operator-(difference_type n) const
		{
			difference_type reindex = index_ - n;
			check_range(reindex);
			return iterator(*outclass_, reindex);
		}
		difference_type operator-(const iterator& other) const
		{
			return index_ - other.index_;
		}
		reference operator[](difference_type n) const
		{
			difference_type reindex = index_ + n;
			check_range(reindex);
			if (reindex == outclass_->size())
				throw std::out_of_range("out of range");
			return (*outclass_)[reindex];
		}
		auto operator<=>(const iterator& other) const { return (*data_)[index_] <=> (*other.data_)[index_]; }
		bool operator==(const iterator& other) const { return (*outclass_)[index_] == (*other.outclass_)[other.index_]; }
		iterator& operator=(const iterator& other)
		{
			index_ = other.index_;
			outclass_ = other.outclass_;
			return *this;
		}
		iterator& operator=(iterator&& other) noexcept
		{
			index_ = other.index_;
			outclass_ = std::move(other.outclass_);
			other.index_ = 0;
			other.outclass_ = nullptr;
			return *this;
		}
	};
	iterator begin() { return iterator(*this, 0); }
	iterator end() { return iterator(*this, size()); }
	VectoCirc(VectoCirc&) noexcept;
	VectoCirc(VectoCirc&&) noexcept;
	VectoCirc(std::vector<T>& WAVE);
	VectoCirc(std::vector<T>&& WAVE);
	VectoCirc() : data_(), start_(0), end_(0){}
	VectoCirc& operator=(VectoCirc&&) noexcept;
	T& operator[](size_t index);
	const T& operator[](size_t index) const;
	operator std::vector<T>()
	{
		std::vector<T> right_order(size());
		copy(begin(), end(), right_order.begin());
		return right_order;
	}
	int frontmove() const;
	int endmove() const;
	void Strelief() const;
	void extend_front(size_t n);
	void extend_back(size_t n);
	void shrink_front(size_t n);
	void shrink_back(size_t n);
	void push_front(const VectoCirc<T>& array);
	void push_front(const std::vector<T>& array);
	void push_front(const T& element);
	void push_back(const VectoCirc<T>& array);
	void push_back(const std::vector<T>& array);
	void push_back(const T& element);
	void slide_front(const std::vector<T>& array);
	void slide_back(const std::vector<T>& array);
	T pop_front();
	T pop_back();
	size_t size() const;
	void clear();
};

template <typename T>
class TailVectoCirc
{
private:
	bool virtualflag = false;
	size_t reservetail = 0, virtualfront = 0;
	VectoCirc<T> vectorcirc;
public:
	class iterator
	{
	private:
		long long index_;
		TailVectoCirc<T>* outclass_;
		void check_range(long long idex) const
		{
			if (idex < 0 || idex > outclass_->vectorcirc.size())
				throw std::out_of_range("out of range");
		}
		constexpr long long real_index(long long id) const
		{
			return id - outclass_->virtualfront;
		}
	public:
		using iterator_category = std::random_access_iterator_tag;
		using value_type = T;
		using difference_type = std::ptrdiff_t;
		using pointer = T*;
		using reference = T&;
		iterator(TailVectoCirc<T>& parent, size_t index) : outclass_(&parent), index_(index) {}
		iterator(const iterator& other) : outclass_(other.outclass_), index_(other.index_) {}
		iterator(iterator&& other) noexcept : outclass_(std::move(other.outclass_)), index_(other.index_)
		{
			other.index_ = 0;
			other.outclass_ = nullptr;
		}
		reference operator*()
		{
			return outclass_->vectorcirc[real_index(index_)];
		}
		iterator& operator++()
		{
			++index_;
			check_range(real_index(index_));
			return *this;
		}
		iterator operator++(int)
		{
			iterator tmp = *this;
			++index_;
			check_range(real_index(index_));
			return tmp;
		}
		iterator& operator--()
		{
			--index_;
			check_range(real_index(index_));
			return *this;
		}
		iterator operator--(int)
		{
			iterator tmp = *this;
			--index_;
			check_range(real_index(index_));
			return tmp;
		}
		iterator operator+(difference_type n) const
		{
			difference_type reindex = index_ + n;
			check_range(real_index(reindex));
			return iterator(*outclass_, reindex);
		}
		iterator operator-(difference_type n) const
		{
			difference_type reindex = index_ - n;
			check_range(real_index(reindex));
			return iterator(*outclass_, reindex);
		}
		difference_type operator-(const iterator& other) const
		{
			return index_ - other.index_;
		}
		reference operator[](difference_type n) const
		{
			difference_type reindex = index_ + n;
			check_range(real_index(reindex));
			if (reindex - outclass_->virtualfront == outclass_->vectorcirc.size())
				throw std::out_of_range("out of range");
			return (*outclass_)[real_index(reindex)];
		}
		auto operator<=>(const iterator& other) const { return outclass_->vectorcirc[real_index(index_)] <=> other.outclass_->vectorcirc[other.real_index(other.index_)]; }
		bool operator==(const iterator& other) const { return outclass_->vectorcirc[real_index(index_)] == other.outclass_->vectorcirc[other.real_index(other.index_)]; }
		iterator& operator=(const iterator& other)
		{
			index_ = other.index_;
			outclass_ = other.outclass_;
			return *this;
		}
		iterator& operator=(iterator&& other) noexcept
		{
			index_ = other.index_;
			outclass_ = std::move(other.outclass_);
			other.index_ = 0;
			other.outclass_ = nullptr;
			return *this;
		}
	};
	iterator begin() { return iterator(*this, 0); }
	iterator end() { return iterator(*this, size()); }
	TailVectoCirc() :vectorcirc() {}
	TailVectoCirc(TailVectoCirc& other);
	TailVectoCirc(TailVectoCirc&& other) noexcept;
	TailVectoCirc& operator=(TailVectoCirc&& other) noexcept;
	TailVectoCirc(VectoCirc<T>& array, size_t vf = 0);
	TailVectoCirc(VectoCirc<T>&& array, size_t vf = 0);
	operator VectoCirc<T>()
	{
		if (virtualflag)
			throw std::logic_error("virtualflag is true");
		return vectorcirc;
	}
	void virtual_transform(size_t rt);
	void virtual_transform();
	T& operator[](size_t index);
	const T& operator[](size_t index) const;
	int frontmove() const;
	int endmove() const;
	void Strelief() const;
	void extend_back(size_t n);
	void shrink_front(size_t n);
	void shrink_back(size_t n);
	void push_back(const VectoCirc<T>& array);
	void push_back(const std::vector<T>& array);
	void push_back(const T element);
	void slide_back(const std::vector<T>& array);
	T pop_front();
	T pop_back();
	size_t size() const;
	void clear();
};