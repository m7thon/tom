/**
 * @file   Sequence.h
 * @author Michael Thon
 *
 */

#ifndef SEQUENCE_H
#define SEQUENCE_H

namespace tom {

constexpr long NoIndex = std::numeric_limits<long>::max();

SWIGCODE(%ignore SequenceData);
class SequenceData {
	friend class cereal::access;
public:
	SequenceData(int nO = 0, int nU = 0) : nO_(nO), nU_(nU) {}
	SequenceData(const std::vector<Symbol>& seq, int nO = 0, int nU = 0) : nO_(nO), nU_(nU), seq_(seq) {}
	SequenceData(unsigned int size, int nO = 0, int nU = 0) : nO_(nO), nU_(nU), seq_(size, 0) {}
  int nO_;                  ///< The size of the output alphabet
	int nU_;                  ///< The size of the input alphabet, or 0 if there are no inputs
	std::vector<Symbol> seq_; ///< the underlying sequence data
};

SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Sequence::repr;)
SWIGCODE(%feature("python:slot", "sq_length", functype="lenfunc") Sequence::size;)


/**
 * This class represents a sequence to be used with the OOM algorithms.\ It may represent a sequence, subsequence view or io-sequence, and stores information about the output and possibly input alphabet.
 */
class Sequence {
	friend class SequenceData;
public:
	typedef Symbol value_type;

/** @name Constructors */
//@{

	/**
	 * Construct a \a Sequence with output alphabet size \a nO and input alphabet size \a nU from a given \a data vector.\ The \a data vector is copied, and the sequence is viewed as an input-output sequence if \a nU != 0.
	 */
	Sequence(const std::vector<Symbol>& data, int nO, int nU = 0);

	/**
	 * Construct a zero \a Sequence with output alphabet size \a nO and input alphabet size \a nU of a given \a length.\ The size of the \a Sequence will be 2 * \a length if it is an input-output sequence, i.e., if \a nU != 0.
	 */
	Sequence(unsigned long length = 0, int nO = 0, int nU = 0);

  /** Construct a \a Sequence from the string representation given by the \a sequenceStr\. This must correspond to what the \a toString() member function produces. */
  Sequence(const std::string& sequenceStr) { std::stringstream iss(sequenceStr); *this << iss; }
//@}

/** @name Accessors */
//@{
	/** Return the size of the input alphabet.\ If this is zero, then the \a Sequence is an ordinary sequence, else an io-sequence. */
	int nU() const { return data_->nU_; }

	/** Return the size of the output alphabet. */
	int nO() const { return data_->nO_; }

	/** Return the size of the represented (sub)-sequence.\ Note that for an io-sequence \f$ u_0o_0\ldots u_{N-1}o_{N-1}\f$ this is 2*N. */
	unsigned long size() const { return labs(size_); }

	/** Return \c true if this is a reversed \a Sequence. */
	bool isReversed() const { return ( size_ < 0 ); }

  /** Return the index of the first symbol of the represented (sub)-sequence in the underlying data vector. */
	unsigned long pos() const { return pos_; }

	/** Return the (input or output) symbol at position \a n. */
	const Symbol at(unsigned long n) const;

	#ifndef SWIG
	/** Access the (input or output) symbol at position \a n. */
	Symbol& at(unsigned long n);

	/** Return the (input or output) symbol at position \a n. */
	const Symbol operator[] (unsigned long n) const { return at(n); }

	/** Access the (input or output) symbol at position \a n. */
	Symbol& operator[] (unsigned long n) { return at(n); }
	#endif
//@}

/** @name Functionality */
//@{
	/** Return a deep copy of this \a Sequence, i\.e\., the copy will use its own memory. */
	Sequence copy() const;

	/** Return an input-output \a Sequence formed by adding the values provided in the output-only \a Sequence \a inputSequence as inputs to the current output \a Sequence\. If the provided \a inputSequence is shorter than this \a Sequence, the remaining inputs will be set to zero. */
	Sequence mergeToIo(const Sequence& inputSequence = Sequence()) const;

	/** Reverse the \a Sequnece.\ Note that this does not actually change the underlying data. */
	void reverse() { size_ = -size_; }

	/** Return a subsequence starting at the given position \a pos and of the given \a size. */
	Sequence sub(unsigned long pos, unsigned long size, bool reverse = false) const;

	/** Return a slice from the \a begin position up to (and not including) the \a end position, or if \a reverse is set to \c true, backwards from the \a end position up to (and not including) the \a begin position\. This follows python syntax, i.e., negative position values are counted from the end of the \a Sequence, with 0 being the first and -1 being the last position\. Furthermore, \a begin and \a end may be set to \a NoIndex, and then extend to the beginning or end of the \a Sequence. **/
	Sequence slice(long begin = NoIndex, long end = NoIndex, bool reverse = false) const;

	/** Alias for sub(unsigned long pos, unsigned long size) */
	Sequence substr(unsigned long pos, unsigned long size) const { return sub(pos, size); }

	/** Return true if the given \a Sequence \a seq is equal to this \a Sequence */
	bool operator ==(const Sequence& seq) const {
		if ((nU() == 0 and seq.nU() > 0) or (nU() > 0 and seq.nU() == 0) or (size() != seq.size())) return false;
		for (unsigned int i = 0; i < size(); ++i) if (at(i) != seq.at(i)) return false;
		return true;
	}

	SWIGCODE(%ignore operator++;)
	/** Increment the first position, which allows using a \a Sequence as an iterator */
	void operator ++() {
		if (size_ == 0) return;
		if (isReversed()) { size_++; }
		else { pos_++; size_--; }
	}
	
	/** Count the number of occurrences of the given \a Sequence \a seq as a sub-sequence of this \a Sequence. */
	unsigned int count(const Sequence& seq) const {
		if ((nU() == 0 and seq.nU() > 0) or (nU() > 0 and seq.nU() == 0)) return 0;
		if (seq.size() == 0) return length();
		unsigned int c = 0;
		for (unsigned int i = 0; i <= size() - seq.size(); i += nU() == 0 ? 1 : 2) if (seq == (sub(i, seq.size()))) c++;
		return c;
	}

//@}

/** @name Special functions for input-outpus sequences */
//@{
	/** Return \c true if this is an input-output sequence, i.e., if the input alphabet size \a nU is non-zero. */
	bool isStrictlyIO() const { return ( nU() != 0 ); }

  /** Return \c true if this sequence is aligned.\ This is always true for normal (non-io) sequences.\ For io-sequences this means that the sequence begins with an input-output pair, i.e., this \a Sequence is aligned if it is not reversed and begins with an input symbol, or is reversed and begins with an output symbol. */
	bool isAlignedIO() const { return ((nU() == 0) or (pos_ % 2 == 0)); }

	/** Return \c true if this sequence is a valid (io)-sequence.\ This is always true for normal (non-io) sequences.\ For io-sequences this means that the sequence is aligned and has an even number of symbols. */
	bool isValidIO() const { return ((nU() == 0) or ((size_%2==0) and isAlignedIO())); }

	/** Return a subsequence starting at the given index \a posIO and of the given \a length, i.e., if this \a Sequence is \f$ u_0o_0\ldots u_{N-1}o_{N-1}\f$, then return a subsequence view to \f$ u_{posIO}o_{posIO}\ldots u_{posIO+length-1}o_{posIO+length-1}\f$.\ For normal (non-io) sequences this is the same as \a sub. */
	Sequence subIO(unsigned long posIO, unsigned long length) const {	return sub( (nU() == 0 ? posIO : 2*posIO), (nU() == 0 ? length : 2*length) ); }

  /** Return the length of this sequence, which is its size for normal (non-io) sequences and half its size for io-sequences.\ Note that this only makes sense for valid (io)-sequences. */
	unsigned long length() const { return (nU() == 0 ? size() : size()/2); }

	/** Return the \a n-th output symbol.\ Note that this is only correct for aligned (io)-sequences. */
	const Symbol o(unsigned long n) const;
	/** Set the \a n-th output symbol.\ Note that this is only correct for aligned (io)-sequences. */
	void o(unsigned long n, Symbol o);

	/** Return the \a n-th input symbol or 0 if this is an output-only \a Sequence.\ Note that this is only correct for aligned and strictly io-sequences. */
	const Symbol u(unsigned long n) const;
	/** Set the \a n-th input symbol.\ Note that this is only correct for aligned and strictly io-sequences. */
	void u(unsigned long n, Symbol u);
//@}

/** @name IO-functions */ //@{
  /** return a \a std::string representation that can be used for saving to file etc. This function should just call the output stream operator. */
	std::string toString() const { std::stringstream oss; *this >> oss; return oss.str(); }
	INSERT_JSON_IO_FUNCTIONS()
	/** return a representation to display in interactive python. */
	std::string repr() const;
  /** read from the given input stream and initialize\. The format must correspond to what the output function produces. */
private:
#ifndef SWIG
	template<class Archive>
  void save(Archive & ar) const {
		const std::string type = "SEQUENCE";
		ar(cereal::make_nvp("Type", type));
		ar(cereal::make_nvp("size", (std::int64_t)(size_)));
		ar(cereal::make_nvp("nO", nO()));
		ar(cereal::make_nvp("nU", nU()));
		ar.setNextName("data"); ar.startNode(); ar.makeArray(); ar.writeName();
		for (unsigned long i = 0; i < size(); ++i) {
			// Note that we must save the data in the original (non-reversed) order!
			ar.saveValue(data_->seq_.at( pos() + i ));
		}
		ar.finishNode();
  }
  template<class Archive>
  void load(Archive & ar) {
		pos_ = 0;
		std::string type;
		ar(cereal::make_nvp("Type", type));
		int nO, nU;
		bool haveSize = true;
		try { std::int64_t len; ar(cereal::make_nvp("size", len)); size_ = len; }
		catch(...) { haveSize = false; }
		ar(cereal::make_nvp("nO", nO));
		ar(cereal::make_nvp("nU", nU));
		ar.setNextName("data"); ar.startNode();
		if (!haveSize) { cereal::size_type s; ar.loadSize(s); size_ = s; }
		data_ = std::make_shared<SequenceData>(size(), nO, nU);
		for (unsigned long i = 0; i < size(); ++i) { ar.loadValue(data_->seq_.at(i)); }
		ar.finishNode();
	}
#endif // SWIG
	std::istream & operator<<(std::istream &istream);
  /** write to the given output stream. */
	std::ostream& operator>>(std::ostream &os) const;
	/** an output helper function */
	std::ostream& writeFormattedData(std::ostream& ostream, unsigned int maxOut = 0) const;
//@}

private:
	unsigned long pos_;             ///< The index position in the \a data_ corresponding to the beginning of this \a Sequence
	long size_;                     ///< The size (in Symbols) of this \a Sequence.\ Negative values indicate a reversed sequence
	std::shared_ptr<SequenceData> data_; ///< a pointer to the underlying \a SequenceData
}; // class Sequence

// The following is some magic to make Sequence objects slice-able and iterable in python
struct stop_iteration {};
#ifdef SWIG
%extend Sequence {
	%typemap(in) PySliceObject* {
		if (!PySlice_Check($input)) { %argument_fail(SWIG_TypeError, "$type", $symname, $argnum); }
		$1 = (PySliceObject *) $input;
	}
	%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) PySliceObject* { $1 = PySlice_Check($input); }
	%feature("python:slot", "mp_subscript", functype="binaryfunc") __getitem__;
	Sequence __getitem__(PySliceObject *slice) throw (std::invalid_argument) {
		Py_ssize_t start, stop, step, sliceLen;
		PySlice_GetIndicesEx(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &start, &stop, &step, &sliceLen);
		if (step == 1) { return self->sub(start, sliceLen); }
		else if (step == -1) { return self->sub(stop+1, sliceLen, true); }
		else { throw std::invalid_argument("slice step must be -1 or 1 for tom::Sequence objects."); }
	}
	Symbol __getitem__(long i) throw (std::out_of_range) {
		if (i < 0) i = i + self->size();
		if ((i < 0) or (i >= self->size())) { throw std::out_of_range("Index out of bounds"); }
		return self->at(i);
	}
	%feature("python:slot", "mp_ass_subscript", functype="objobjargproc") __setitem__;
	void __setitem__(long i, const Symbol& val) throw (std::out_of_range) {
		if (i < 0) i = i + self->size();
		if ((i < 0) or (i >= self->size())) { throw std::out_of_range("Index out of bounds"); }
		self->at(i) = val;
	}
	%feature("python:slot", "tp_iter", functype="getiterfunc") __iter__;
	Sequence __iter__() { return tom::Sequence(*self); }
	%typemap(throws) stop_iteration {
    (void)$1;
    SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
    SWIG_fail;
  }
	%catches(stop_iteration) __next__();
	%feature("python:slot", "tp_iternext", functype="iternextfunc") __next__;
	Symbol __next__() {
		if (self->size() == 0) { throw tom::stop_iteration(); }
		tom::Symbol ret = self->at(0);
		++(*self);
		return ret;
	}	
};
#endif // SWIG

/** A vector of sequences */
typedef std::vector<Sequence> Sequences;

#ifndef SWIG
// ##################################################################
//                         IMPLEMENTATION
// ##################################################################

//Constructors and Destructors
Sequence::Sequence(const std::vector<Symbol>& data, int nO, int nU) : pos_(0), size_(0) {
	data_ = std::make_shared<SequenceData>(data, nO, nU);
	size_ = data_->seq_.size();
}

Sequence::Sequence(unsigned long length_, int nO, int nU) : pos_(0), size_(0) {
	size_ = (nU == 0 ? length_ : 2 * length_);
	data_ = std::make_shared<SequenceData>(size_, nO, nU);
}

//Accessors
const Symbol Sequence::at(unsigned long n) const {
	assert(n < size());
	return data_->seq_.at( isReversed() ? pos() + size() - 1 - n : pos() + n );
}
Symbol& Sequence::at(unsigned long n) {
	assert(n < size());
	return data_->seq_.at( isReversed() ? pos() + size() - 1 - n : pos() + n );
}

//Utilities
Sequence Sequence::sub(unsigned long pos, unsigned long size, bool reverse) const {
	assert((pos < this->size()) and (pos + size <= this->size()));
	Sequence seq(*this);
	seq.pos_ = ( isReversed() ? this->pos() + this->size() - pos - size : this->pos() + pos );
	seq.size_ = ( isReversed() ? -size : size );
	if (reverse) { seq.reverse(); }
	return seq;
}

Sequence Sequence::slice(long begin, long end, bool reverse) const {
	if (size() == 0) return Sequence();
	Sequence seq(*this);
	if (not reverse) {
		if (begin == NoIndex) { begin = 0; }
		if (begin < 0) { begin = std::max(0L, begin + long(size())); }
		if (end == NoIndex) { end = size(); }
		if (end < 0) { end += long(size()); }
	} else { // reverse:
		if (begin == NoIndex) { begin = long(size()) - 1; }
		if (begin < 0) { begin += size(); }
		if (end < 0) { end = std::max(-1L, long(size()) + end); }
		if (end == NoIndex) { end = -1; }
		std::swap(++begin, ++end);
	}
	seq.pos_ =  ( isReversed() ? long(pos()) + long(size()) - end : long(pos()) + begin );
	seq.size_ = ( isReversed() ? -std::max(0L, end - begin) : std::max(0L, end - begin) );
	if (reverse) { seq.reverse(); }
	return seq;
}


Sequence Sequence::copy() const {
	Sequence seq(size(), nO(), nU());
	for (unsigned long n = 0; n < length(); ++n) {
		seq.o(n, o(n));
		if (nU() != 0) seq.u(n, u(n));
	}
	return seq;
}

Sequence Sequence::mergeToIo(const Sequence& inputSequence) const {
	unsigned long newnU = inputSequence.nO();
	if (newnU == 0) { newnU = 1; }
	Sequence seq(length(), nO(), newnU);
	unsigned long len = length();
	for (unsigned long n = 0; n < len; ++n) { seq.o(n, o(n)); }
	len = std::min(len, inputSequence.length());
	for (unsigned long n = 0; n < len; ++n) { seq.u(n, inputSequence.o(n)); }
	return seq;
}


//Special IO-Sequence functions
const Symbol Sequence::o(unsigned long n) const {
	if (nU() == 0) { return at(n); }
	else { return ( isReversed() ? at(2*n) : at(2*n+1) ); }
}
void Sequence::o(unsigned long n, Symbol o) {
	if (nU() == 0) { at(n) = o; }
	else { if (isReversed()) at(2*n) = o; else at(2*n+1) = o; }
}
const Symbol Sequence::u(unsigned long n) const {
	if (nU() == 0) return 0;
	return ( isReversed() ? at(2*n+1) : at(2*n) );
}
void Sequence::u(unsigned long n, Symbol u) {
	if (nU() == 0) return;
	if (isReversed()) at(2*n+1) = u; else at(2*n) = u;
}


//IO-functions
std::istream & Sequence::operator<<(std::istream &istream) {
	pos_ = 0;
	std::string description;
	istream >> description >> description;
	istream >> size_;
	int nU, nO;
	istream >> description >> nU;
	istream >> description >> nO;
	if (nU != 0) size_ *= 2;
	data_ = std::make_shared<SequenceData>(size_, nO, nU);
	for (unsigned long n = 0; n < size_; ++n)
		istream >> data_->seq_.at(n);
	return istream;
}

std::ostream& Sequence::writeFormattedData(std::ostream& ostream, unsigned int maxOut) const {
	if ((maxOut == 0) or (maxOut > size())) maxOut = size();
	bool align = isAlignedIO();
	if (maxOut != 0) {
		if (nU() != 0) {
			for (unsigned long n = 0; n < maxOut-1; ++n) {
				align = not align;
				ostream << at(n) << ( align ? "  " : " " );
			}
		}
		else
			for (unsigned long n = 0; n < maxOut-1; ++n)
				ostream << at(n) << " ";
		ostream << at(maxOut-1);
	}
	return ostream;
}

std::ostream& Sequence::operator>>(std::ostream &ostream) const {
  ostream << "SEQUENCE " << "Length: " << length() << " nU: " << nU() << " nO: " << nO() << std::endl;
	writeFormattedData(ostream);
	return ostream;
}

std::string Sequence::repr() const {
	std::stringstream os;
	os << "SEQUENCE(" << length() << "," << nO() << "," << nU() << "): ";
	writeFormattedData(os, nU() != 0 ? 20 : 10);
	if (length() > 10) os << " ...";
	return os.str();
}

#endif // SWIG

} // tom

#endif // SEQUENCE_H
