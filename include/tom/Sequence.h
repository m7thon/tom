/**
 * @file   Sequence.h
 * @author Michael Thon
 *
 */

#ifndef SEQUENCE_H
#define SEQUENCE_H

namespace tom {

SWIGCODE(%feature("docstring") Sequence
"Sequence(length = 0, nO = 0, nU = 0)\n"
"Sequence(seq, nO, nU = 0)\n"
"Sequence(json_representation)\n"
"\n"
"This object represents a sequence, subsequence view or io-sequence and stores the\n"
"size `nO` of the output and `nU` of the input alphabet. If the size of the input\n"
"alphabet is zero, this is just an ordinary sequence o_0...o_{N-1} of symbols o_t\n"
"with zero-based indexing.\n"
"\n"
"An io-sequence is represented as a simple sequence u_0o_0...u_{N-1}o_{N-1} of\n"
"inputs u_t and outputs o_t. For io-sequences, we distinguish `size` and `length`:\n"
"For the io-sequence u_0o_0...u_{N-1}o_{N-1}, the `length` is N, while the `size`\n"
"is 2N, and the symbol `at(n)` is `u(n/2)` if n is even and `o((n-1)/2)` if n is\n"
"odd, or, conversely, the n-th input symbol `u(n)` is the symbol `at(2n)`, and\n"
"the n-th output symbol `o(n)` is the symbol `at(2n+1)`. For standard sequences,\n"
"`o(n)` is always `at(n)` and `u(n)` is always zero.\n"
"\n"
"This object always represents a view to underlying sequence data, i.e., copies,\n"
"slices and subsequences always point to the same underlying data. To obtain a\n"
"real deep copy, the `copy` member function is provided.\n"
"\n"
"Constructors\n"
"------------\n"
"Sequence(length = 0, nO = 0, nU = 0) -> a Sequence of given `length`, output\n"
"    alphabet size `nO` and input alphabet size `nU` initialized with zeros. In\n"
"    the case of an io-sequence (if `nU` != 0), the `size` will be 2 * `length`.\n"
"Sequence(seq, nO, nU = 0) -> a Sequence of given output alphabet size `nO` and\n"
"    input alphabet size `nU` constructed from the given `seq`, which may be a\n"
"    list [u_0, ..., u_{N-1}] (or [u_0, o_0, ..., u_{N-1}, o_{N-1}] for an io-\n"
"    sequence), or a std::vector<int>. The contents of `seq` will be copied.\n"
"Sequence(json_representation) -> a Sequence corresponding to the given string\n"
"    `json_representation`. The format should correspond to what `.toJSON()`\n"
"    produces.\n"
"\n";)

constexpr long NoIndex = std::numeric_limits<long>::min();

SWIGCODE(%ignore SequenceData);
class SequenceData {
	friend class cereal::access;
	friend class Sequence;
public:
	SequenceData(Symbol nO = 0, Symbol nU = 0) : nO_(nO), nU_(nU) {}
	SequenceData(const std::vector<Symbol>& seq, Symbol nO = 0, Symbol nU = 0) : nO_(nO), nU_(nU), seq_(seq) {}
	SequenceData(unsigned long size, Symbol nO = 0, Symbol nU = 0) : nO_(nO), nU_(nU), seq_(size, 0) {}
    Symbol nO_;                  ///< The size of the output alphabet
	Symbol nU_;                  ///< The size of the input alphabet, or 0 if there are no inputs
	std::vector<Symbol> seq_;    ///< the underlying sequence data
};

SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Sequence::repr;)
SWIGCODE(%feature("python:slot", "sq_length", functype="lenfunc") Sequence::size;)
/**
 * This is the basic class to represent a sequence, subsequence view or io-sequence, and stores information about the size of the output and input alphabet\. If the size of the input alphabet is zero, this is just an ordinary sequence of symbols \f$ o_0\ldots o_{N-1}\f$\, with zero-based indexing\. Note that an io-sequence is represented as a simple sequence of inputs \f$u_t\f$ and outputs \f$o_t\f$, i.e., as \f$ u_0o_0\ldots u_{N-1}o_{N-1}\f$\. For such an io-sequence, the \a length is \a N, while the \a size is \a 2N\. The symbol \a at(n) is \f$u_{n/2}\f$ if \a n is even and \f$o_{(n-1)/2}\f$ if \a n is odd\. Conversely, the \a n-th input symbol \a u(n) is the symbol \a at(2n), and the \a n-th output symbol \a o(n) is the symbol \a at(2n+1).
 */
class Sequence {
public:
	typedef Symbol value_type;

/** @name Constructors */
//@{
    /** Construct a \a Sequence with output alphabet size \a nO and input alphabet size \a nU from a given \a data vector.\ The \a data vector is copied, and the sequence is viewed as an input-output sequence if \a nU != 0. */
	Sequence(const std::vector<Symbol>& data, Symbol nO, Symbol nU = 0);

    /** Construct a \a Sequence of zeros with output alphabet size \a nO and input alphabet size \a nU of a given \a length.\ The size of the \a Sequence will be 2 * \a length if it is an input-output sequence, i.e., if \a nU != 0. */
    Sequence(unsigned long length = 0, Symbol nO = 0, Symbol nU = 0);

    /** Construct a \a Sequence from the given \a json_representation\. This must correspond to what the \a toJSON() member function produces. */
    Sequence(const std::string& json_representation) { assert(false, "Not implemented"); }
//@}
    
/** @name Accessors and Properties */
//@{
    /** Return the size of the input alphabet.\ If this is zero, then the \a Sequence is an ordinary sequence, else an io-sequence. */
	Symbol nU() const { return data_->nU_; }
    
    /** Return the size of the output alphabet. */
	Symbol nO() const { return data_->nO_; }
    
    /** Return the size of the represented sequence.\ Note that for an io-sequence this is the sum of the number of input symbols and output symbols.\ Generally, use \a length() instead. */
	unsigned long rawSize() const { return labs(size_); }

    /** Return the length of this sequence\. For io-sequences this is the number of input-output symbols (i.e., time steps), which is generally half the size. */
	unsigned long length() const { return nU() == 0 ? rawSize() : isFrontAligned() ? (rawSize() + 1) / 2 : rawSize() / 2 + 1; }
    
    /** Return the output symbol at the \a n-th position. */
    Symbol o(long n) const { return data_->at(indexToPos(n, true)); }

    /** Set the output symbol at the \a n-th posiiton. */
	void o(long n, Symbol o) { data_->at(indexToPos(n, true)) = o; }

    /** Return the input symbol at the \a n-th position, or zero if this is an output-only sequence. */
    Symbol u(long n) const { if (nU() == 0) return 0; return data_->at(indexToPos(n, false)); }
    
    /** Set the output symbol at the \a n-th posiiton, or do nothing if this is an output-only sequence. */
	void u(long n, Symbol u) { if (nU() == 0) return; data_->at(indexToPos(n, false)) = u; }

    /** Return \c true if this is a reversed \a Sequence\. This is relevant for io-sequences, since the order of input and outputs is also reversed.*/
	bool isReversed() const { return size_ < 0; }
	
    bool isFrontAligned() const { return nU() == 0 or pos_ % 2 == 0; }
	bool isBackAligned() const { return nU() == 0 or (pos_ + rawSize()) % 2 == 0; }
	bool isAligned() const { return isFrontAligned() and isBackAligned(); }
//@}
    
/* Functionality */
//@{
	/** Return a deep copy of this \a Sequence, i\.e\., the copy will use its own memory. */
	Sequence copy() const;

	/** Return an input-output \a Sequence formed by adding the values provided in the output-only \a Sequence \a inputSequence as inputs to the current output \a Sequence\. If the provided \a inputSequence is shorter than this \a Sequence, the remaining inputs will be set to zero. */
	Sequence mergeToIo(const Sequence& inputSequence = Sequence()) const;

	/** Reverse the \a Sequnece.\ Note that this does not actually change the underlying data. */
	void reverse() { size_ = -size_; }

	/** Return a subsequence starting at the given position \a pos and of the given \a size. */
	Sequence sub(unsigned long pos, unsigned long size, bool reverse = false) const __attribute__((deprecated));

	/** Return a slice from the \a begin position up to (and not including) the \a end position, or if \a reverse is set to \c true, backwards from the \a end position up to (and not including) the \a begin position\. This follows python syntax, i.e., negative position values are counted from the end of the \a Sequence, with 0 being the first and -1 being the last position\. Furthermore, \a begin and \a end may be set to \a NoIndex, and then extend to the beginning or end of the \a Sequence. **/
	Sequence slice(long begin = NoIndex, long end = NoIndex, bool reverse = false) const __attribute__((deprecated));

	/** Return true if the given \a Sequence \a seq is equal to this \a Sequence */
	bool operator ==(const Sequence& seq) const {
		if ((nU() == 0 and seq.nU() > 0) or (nU() > 0 and seq.nU() == 0) or (rawSize() != seq.rawSize()) or (isFrontAligned() != seq.isFrontAligned())) return false;
		if (nU() == 0 or ((isReversed() == seq.isReversed()))) {
			for (unsigned int i = 0; i < rawSize(); ++i) if (rawAt(i) != seq.rawAt(i)) return false;
			return true;
		} /* else { // can a reversed io-sequence equal a non-reversed one?
			if (isAligned()) {
				for (unsigned int i = 0; i < length(); ++i) if ((o(i) != seq.o(i)) or u(i) != seq.u(i)) return false;
				return true;
			}
		} */
		return false;
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
        for (unsigned int i = 0; i <= length() - seq.length(); ++i) {
            if (seq == (sub(i, seq.length()))) c++;
        }
		return c;
	}
//@}

/** @name Special functions for input-outpus sequences */
//@{
	/** Return \c true if this is an input-output sequence, i.e., if the input alphabet size \a nU is non-zero. */
	bool isIO() const { return ( nU() != 0 ); }

	/** Return a subsequence starting at the given index \a posIO and of the given \a length, i.e., if this \a Sequence is \f$ u_0o_0\ldots u_{N-1}o_{N-1}\f$, then return a subsequence view to \f$ u_{posIO}o_{posIO}\ldots u_{posIO+length-1}o_{posIO+length-1}\f$.\ For normal (non-io) sequences this is the same as \a sub. */
	Sequence subIO(unsigned long posIO, unsigned long length) const {	return sub( (nU() == 0 ? posIO : 2*posIO), (nU() == 0 ? length : 2*length) ); }

  /** Return the length of this sequence, which is its size for normal (non-io) sequences and half its size for io-sequences.\ Note that this only makes sense for valid (io)-sequences. */
	unsigned long length() const { return (nU() == 0 ? size() : size()/2); }
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
	long indexToPos(long n = 0, bool forOutputSymbol = true) const {
		long p = forOutputSymbol;
		if (n<0) n += length();
		if (nU() == 0) { p += isReversed() ? pos_ + rawSize() - 1 - n : pos_ + n; }
		else { p += isReversed() ? 2*((pos_ + rawSize() - 1)/2 - n) : 2*(pos_/2 + n); }
		if (p < pos_ or p >= pos_ + rawSize()) throw std::out_of_range();
		return p;
	}
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
	%feature("python:slot", "nb_nonzero", functype="inquiry") __nonzero__;
	bool __nonzero__() const { return (self->nO() != 0); }
	bool __bool__() const { return (self->nO() != 0); }
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
	%typemap(throws) tom::stop_iteration {
    (void)$1;
    SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
    SWIG_fail;
  }
	%feature("python:slot", "tp_iternext", functype="iternextfunc") __next__;
	Symbol __next__() throw(tom::stop_iteration) {
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
	Sequence seq(*this);
	if (not reverse) {
		if (begin == NoIndex) { begin = 0; }
		else {
			if (begin < 0) { begin = std::max(0L, begin + long(size())); }
			if (begin > long(size())) { begin = long(size()); }
		}
		if (  end == NoIndex) {   end = long(size()); }
		else {
			if (  end < 0) { end   = std::max(0L,   end + long(size())); }
			if (  end > long(size())) {   end = long(size()); }
		}
	} else { // reverse:
		if (begin == NoIndex) { begin = long(size()) - 1; }
		else {
			if (begin < 0) { begin = std::max(-1L, begin + long(size())); }
			if (begin > long(size())-1) { begin = long(size())-1; }
		}
		if (  end == NoIndex) {   end = -1L; }
		else {
			if (  end < 0) {   end = std::max(-1L,   end + long(size())); }
			if (  end > long(size())-1) {   end = long(size())-1; }
		}
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
	os << "Sequence(" << length() << "," << nO() << "," << nU() << "): ";	
	writeFormattedData(os, 10);
	return os.str();
}

#endif // SWIG

} // tom

#endif // SEQUENCE_H
