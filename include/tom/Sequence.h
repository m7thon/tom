#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "tom.h"

namespace tom {
    
constexpr long NoIndex = std::numeric_limits<long>::min();

SWIGCODE(%ignore SequenceData);
class SequenceData {
	friend class cereal::access;
	friend class Sequence;
public:
	SequenceData(Symbol nO = 0, Symbol nU = 0) : nO_(nO), nU_(nU) {}
	SequenceData(const std::vector<Symbol>& seq, Symbol nO = 0, Symbol nU = 0) : nO_(nO), nU_(nU), seq_(seq) {}
	SequenceData(long size, Symbol nO = 0, Symbol nU = 0) : nO_(nO), nU_(nU), seq_(size, 0) {}
    Symbol nO_;                  ///< The size of the output alphabet
	Symbol nU_;                  ///< The size of the input alphabet, or 0 if there are no inputs
	std::vector<Symbol> seq_;    ///< the underlying sequence data
};

SWIGCODE(%feature("python:slot", "tp_repr", functype="reprfunc") Sequence::repr;)
SWIGCODE(%feature("python:slot", "sq_length", functype="lenfunc") Sequence::rawSize;)
/** This object represents a sequence, subsequence view or io-sequence and stores the size `nO()` of the output and `nU()` of the input alphabet. If the size of the input alphabet is zero, this is just an ordinary sequence \f$o_0...o_{N-1}\f$ of symbols \f$o_t\f$ with zero-based indexing. An io-sequence is represented as a simple sequence \f$u_0o_0...u_{N-1}o_{N-1}\f$ of input symbols \f$u_t\f$ and output symbols \f$o_t\f$. For io-sequences, we distinguish *size* and *length*: *size* is always the number of symbols, while the *length* is the number of io symbol pairs, which is just the *size* for ordinary sequences, and 2 * *size* for (aligned) io-sequences.
 
    There are three ways to interact with this sequence:
    1. The `raw...()` methods. These just access the sequence as a raw symbol sequence, treating each input or output symbol as a separate symbol. I.e., for the io sequence \f$u_0o_0...u_{N-1}o_{N-1}\f$, the `rawSize()` is \f$2N\f$, and `rawAt(i)` is \f$u_{i/2}\f$ if \f$i\f$ is even or \f$o_{i/2}\f$ if \f$i\f$ is odd. Etc.
    2. The methods not prefixed by "raw" treat each io-symbol *pair* as one symbol: `at(i)` is the symbol-pair \f$(u_i, o_i)\f$ (actually, the subsequence at index \f$i\f$ of length 1: `sub(i,1)`), the methods `u(i)` and `o(i)` return the input symbol \f$u_i\f$ or respectively output symbol \f$o_i\f$, and the `length()` is the *length* of the sequence -- the number N of io symbol pairs. For ordinary sequences these methods are equivalent to the `raw...()` ones.
    3. Access using python `[]`-syntax (including slicing) and python iteration just treats all sequences as raw symbol sequences (as in 1.).

    This object always represents a view to underlying sequence data, i.e., copies, slices and subsequences always point to the same underlying data. To obtain a real deep copy, the `copy()` member function is provided.
 */
class Sequence {
public:
	typedef Symbol value_type;
    
/** @name Constructors */
//@{
    /** Construct a `Sequence` of given output alphabet size `nO` and input alphabet size `nU` from the given `symbol_list`. This may be a list `[u_0, ..., u_{N-1}]` (or `[u_0, o_0, ..., u_{N-1}, o_{N-1}]` for an io-sequence), or a `std::vector<int>`. The contents of `symbol_list` is copied.
     */
    Sequence(const std::vector<Symbol>& symbol_list, Symbol nO, Symbol nU = 0) {
        data_ = std::make_shared<SequenceData>(symbol_list, nO, nU);
        size_ = data_->seq_.size();
    }

    /** Construct a `Sequence` of given `length`, output alphabet size `nO` and input alphabet size `nU` initialized with zeros. In the case of an io-sequence (if `nU != 0`), the `rawSize()` will be 2 * `length`.
     */
    Sequence(long length = 0, Symbol nO = 0, Symbol nU = 0) {
        size_ = (1 + (nU != 0)) * length;
        data_ = std::make_shared<SequenceData>(rawSize(), nO, nU);
    }

    /** Construct a `Sequence` corresponding to the given string `json_representation`. The format must correspond to what `toJSON()` produces.
     */
    Sequence(const std::string& json_representation) { fromJSON(json_representation.c_str()); }
//@}
    
/** @name Accessors and Properties */
//@{
    /** Return the size of the input alphabet. If this is zero, then this is an ordinary sequence, else an io-sequence.
     */
	Symbol nU() const { return data_->nU_; }
    
    /** Return the size of the output alphabet.
     */
	Symbol nO() const { return data_->nO_; }
    
    /** Return \c true if this is an input-output sequence, i.e., if the input alphabet size \c nU() is non-zero.
     */
    bool isIO() const { return ( nU() != 0 ); }
//@}
    
/** @name Interface as raw symbol sequence */
//@{
    /** Return the *size* of the represented sequence, i.e., the raw symbol count, counting each input and output as one symbol. Generally, use `length()` instead.
     */
	long rawSize() const { return labs(size_); }
    
    /** Return the symbol at index `idx`, treating io-sequences as raw sequences.
     */
    Symbol rawAt(long idx) const CHECK(throw (std::out_of_range)) { return data_->seq_[ checkPos( rawIndexToPos( normalizeRawIndex( idx ) ) ) ]; }
    
    /** Set the symbol at index `idx` to `x`, treating io-sequences as raw sequences. Negative indexing is supported.
     */
    void rawAt(long idx, Symbol x) CHECK(throw (std::out_of_range)) { data_->seq_[ checkPos( rawIndexToPos( normalizeRawIndex( idx ) ) ) ] = x; }

    /** Return a subsequence starting at the given position index `idx` and of the given `size`, treating io-sequences as raw sequences. Negative indexing is supported, and if `size` is negative, a reverse sequence starting at the `idx` is returned.
     */
    Sequence rawSub(long idx, long size) const CHECK(throw (std::out_of_range)) {
        if (size == 0) { Sequence seq(*this); seq.size_ = 0; return seq; }
        idx = normalizeRawIndex( idx );
        long p1 = checkPos( rawIndexToPos( idx ) );
        long p2 = checkPos( rawIndexToPos( idx + size + (size > 0 ? -1 : 1) ) );
        Sequence seq(*this);
        seq.pos_ = std::min(p1, p2);
        seq.size_ = isReversed() ? -size : size;
        return seq;
    }

    /** Return a subsequence from the \c begin index up to (and not including) the \c end index, or if \c forwards is set to \c false, a reverse sequence from the \c begin position up to (and not including) the \c end position. Negative indexing is supported, and \c begin and \c end may be set to \c NoIndex, and then extend to the beginning or end of the sequence, depending on `reverse`.
     */
    Sequence rawSlice(long begin = NoIndex, long end = NoIndex, bool forwards = true) const CHECK(throw (std::out_of_range)) {
        if (begin == NoIndex) { begin = forwards ? 0 : rawSize() - 1; }
        else { begin = normalizeRawIndex( begin ); }
        if (end == NoIndex) { end = forwards ? rawSize() : -1; }
        else {
            end = normalizeRawIndex( end );
            CHECK( if (end < 0) throw std::out_of_range("sequence index out of range"); )
        }
        long size = end - begin;
        CHECK( if ( (begin < 0 and size != 0) or (size < 0 and forwards) or (size > 0 and not forwards) ) { throw std::out_of_range("sequence index out of range"); } )
        return rawSub(begin, size);
    }
//}
    
/** @name Standard interface */
//@{
    /** Return \c true if this is a reversed sequence. This is relevant for io-sequences, since the order of input and outputs is then also reversed.
     */
    bool isReversed() const { return size_ < 0; }
        
    /** Return `true` if this sequence is io-aligned with respect to its beginning in the underlying data, i.e., if either:
        * this is a plain (non-io) sequence
        * this io-sequence is not reversed and begins with an input symbol
        * this io-sequence is reversed and ends with an input symbol
     */
    bool isFrontAligned() const { return nU() == 0 or size_ == 0 or pos_ % 2 == 0; }

    /** Return `true` if this sequence is io-aligned with respect to its end in the underlying data, i.e., if either:
        * this is a plain (non-io) sequence
        * this io-sequence is not reversed and ends with an output symbol
        * this io-sequence is reversed and begins with an output symbol
     */
    bool isBackAligned() const { return nU() == 0 or size_ == 0 or (pos_ + rawSize()) % 2 == 0; }

    /** Return true if this sequence is io-aligned with respect to its underlying data, i.e., if it is front and back aligned.
     */
    bool isAligned() const { return isFrontAligned() and isBackAligned(); }

    /** Return the io symbol pair at index `idx`, where each index covers one io-pair. This returns `sub(idx, 1)`, so even for plain sequences, the return value is not a symbol. For plain sequences, generally use `rawAt(idx)` or `o(idx)` instead.
     
        Negative indexing is supported.
     */
    Sequence at(long idx) const CHECK(throw (std::out_of_range)) { return sub(idx, 1); }
    
    /** Return the output symbol at indx `idx`, where each index covers one io-pair. Negative indexing is supported.
     */
    Symbol o(long idx) const CHECK(throw (std::out_of_range)) { return data_->seq_[ checkPos( indexToPos( normalizeIndex( idx ), true ) ) ]; }

    /** Set the output symbol at index `idx` to `o`, where each index covers one io-pair. Negative indexing is supported.
     */
	void o(long idx, Symbol o) CHECK(throw (std::out_of_range)) { data_->seq_[ checkPos( indexToPos( normalizeIndex( idx ), true ) ) ] = o; }

    /** Return the input symbol at index `idx`, where each index covers one io-pair, or return zero if this is a plain (non-io) sequence. Negative indexing is supported.
     */
    Symbol u(long idx) const CHECK(throw (std::out_of_range)) { if (nU() == 0) return 0; return data_->seq_[ checkPos( indexToPos( normalizeIndex( idx ), false ) ) ]; }
    
    /** Set the input symbol at index `idx` to `u`, where each index covers one io-pair, or do nothing if this is a plain (non-io) sequence. Negative indexing is supported.
     */
	void u(long idx, Symbol u) CHECK(throw (std::out_of_range)) { if (nU() == 0) return; data_->seq_[ checkPos( indexToPos( normalizeIndex( idx ), false ) ) ] = u; }
    
    /** Return the *length* of this sequence. For plain sequences this is the same as `rawSize()`. For (aligned) io-sequences this is the number of io-symbol pairs, which is half the *size*.
     
        For unaligned io-sequences, the *length* means the number of covered io-sequence-pair indices. Example:\n
        For the io-sequence \f$o_0 u_1o_1 ... u_{N-2}o_{N-2} u_{N-1}\f$, which is neither front nor back aligned, the *length* is \f$N\f$, since \f$N\f$ io-pair indices are covered, but the *size* -- the number of raw symbols -- is only \f$2N - 2\f$.
     */
    long length() const { return nU() == 0 ? rawSize() : isFrontAligned() ? (rawSize() + 1) / 2 : rawSize() / 2 + 1; }

    /** Return a subsequence starting at the given position index `idx` and of the given `length`, where each index covers one io-pair. Negative indexing is supported, and if `length` is negative, a reverse sequence starting at the `idx` is returned.
     */
    Sequence sub(long idx, long length) const CHECK(throw (std::out_of_range)) {
        if (nU() == 0) return rawSub(idx, length);
        if (length == 0) { Sequence seq(*this); seq.size_ = 0; return seq; }
        // convert negatives and check:
        idx = normalizeIndex( idx );
        long last_idx;
        bool reverse = length < 0;
        if (reverse) { last_idx = idx; idx = idx + length + 1; length = -length; }
        else { last_idx = idx + length - 1; }
#ifdef TOM_CHECK
        {
            long this_length = this->length();
            if (idx < 0 or idx >= this_length) { throw std::out_of_range("sequence index out of range"); }
            if (last_idx < 0 or last_idx >= this_length) { throw std::out_of_range("sequence index out of range"); }
        }
#endif
        // convert to positions:
        long p1 = indexToPos(      idx,  isReversed() );
        long p2 = indexToPos( last_idx, !isReversed() );
        if (p2 < p1) {
            std::swap(p1, p2);
            reverse = !reverse;
        }
        p1 = std::max(pos_, p1);
        p2 = std::min(pos_ + rawSize() -1, p2);
        Sequence seq(*this);
        seq.pos_ = p1;
        seq.size_ = p2 - p1 + 1;
        if (reverse) { seq.reverse(); }
        return seq;
    }

    /** Return a subsequence from the \c begin index up to (and not including) the \c end index, or if \c forwards is set to \c false, a reverse sequence from the \c begin position up to (and not including) the \c end position, where each index covers one io-pair. Negative indexing is supported, and \c begin and \c end may be set to \c NoIndex, and then extend to the beginning or end of the sequence, depending on `reverse`.
     */
    Sequence slice(long begin, long end = NoIndex, bool forwards = true) const CHECK(throw (std::out_of_range)) {
        if (begin == NoIndex) { begin = forwards ? 0 : length() - 1; }
        else { begin = normalizeIndex( begin ); }
        if (  end == NoIndex) {   end = forwards ? length() : -1; }
        else { end = normalizeIndex( end ); }
        long length = end - begin;
        CHECK( if ( (begin < 0 and length != 0) or (length < 0 and forwards) or (length > 0 and not forwards) ) { throw std::out_of_range("sequence index out of range"); } )
        return sub(begin, length);
    }
//@}
    
/** @name Functionality */
//@{
	/** Return a deep copy of this \c Sequence, i\.e\., the copy will use its own memory.
     */
    Sequence copy() const {
        Sequence seq(length(), nO(), nU());
        if (!isFrontAligned()) { seq.pos_++; }
        seq.size_ = size_;
        for (unsigned long i = 0; i < rawSize(); ++i) { seq.rawAt(i, rawAt(i)); }
        return seq;
    }

    /** Reverse this \c Sequnece.\ Note that this does not actually change the underlying data.
     */
	void reverse() { size_ = -size_; }

	/** Return \c true if the given \c Sequence \c seq is equal to this \c Sequence
     */
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

	SWIGCODE(%ignore incr_as_python_iterator_only;)
	/** Increment the first position, which allows using a \c Sequence as an iterator */
	void incr_as_python_iterator_only() {
        if (isReversed()) { size_ += 1; }
        else { pos_ += 1; size_ -= 1; }
	}
	
	/** Count the number of occurrences of the given \c Sequence \c seq as a sub-sequence of this \c Sequence.
     */
	unsigned int count(const Sequence& seq) const {
		if ((nU() == 0 and seq.nU() > 0) or (nU() > 0 and seq.nU() == 0)) return 0;
		if (seq.length() == 0) return length();
		unsigned int c = 0;
        for (unsigned int i = 0; i <= length() - seq.length(); ++i) {
            if (seq == (sub(i, seq.length()))) c++;
        }
		return c;
	}
//@}

/** @name IO-functions */ //@{
	INSERT_JSON_IO_FUNCTIONS()
	/** return a string representation to display in python. */
	std::string repr() const;
private:
    template<class Archive>
    void save(Archive & ar) const {
        const std::string type = "SEQUENCE";
        ar(cereal::make_nvp("Type", type));
        ar(cereal::make_nvp("size", (std::int64_t)(size_)));
        ar(cereal::make_nvp("nO", nO()));
        ar(cereal::make_nvp("nU", isFrontAligned() ? nU() : -nU()));
        ar.setNextName("data"); ar.startNode(); ar.makeArray(); ar.writeName();
        for (unsigned long i = 0; i < rawSize(); ++i) {
            // Note that we must save the data in the original (non-reversed) order!
            ar.saveValue(data_->seq_.at(pos_ + i ));
        }
        ar.finishNode();
    }
    template<class Archive>
    void load(Archive & ar) {
        std::string type;
        ar(cereal::make_nvp("Type", type));
        int nO = 0;
        int nU = 0;
        size_ = 0;
        bool haveSize = true;
        try { std::int64_t len; ar(cereal::make_nvp("size", len)); size_ = len; }
        catch(...) { haveSize = false; }
        ar(cereal::make_nvp("nO", nO));
        ar(cereal::make_nvp("nU", nU));
        if (nU < 0) { pos_ = 1; nU = -nU; } // the saved sequence is not aligned
        else { pos_ = 0; }
        ar.setNextName("data"); ar.startNode();
        if (!haveSize) { cereal::size_type s; ar.loadSize(s); size_ = s; }
        int data_size = pos_ + rawSize();
        if ((nU != 0) and (data_size % 2 != 0)) { data_size++; }
        data_ = std::make_shared<SequenceData>(data_size, nO, nU);
        for (unsigned long i = 0; i < rawSize(); ++i) { ar.loadValue(data_->seq_.at( pos_ + i )); }
        ar.finishNode();
    }
    std::ostream& writeFormattedData(std::ostream& ostream, long maxOut = 0) const;
//@}

private:
	long pos_ = 0;         ///< The index position in the \c data_ corresponding to the beginning of this \c Sequence
	long size_ = 0;                 ///< The size (in Symbols) of this \c Sequence.\ Negative values indicate a reversed sequence
	std::shared_ptr<SequenceData> data_; ///< a pointer to the underlying \c SequenceData

    /** Return the corresponding raw index for a raw negative index (i.e., treating this as a non-io sequences). */
    long normalizeRawIndex(long idx) const { if (idx < 0) return idx + rawSize(); return idx; }
    /** Return the corresponding io-index for a negative io-index. */
    long normalizeIndex(long idx) const { if (idx < 0) return idx + length(); return idx; }
        
    /** Helper function that returns the position in the underlying data for a given raw normalized index (i.e., treating this as a non-io sequences). */
    long rawIndexToPos(long idx = 0) const { return isReversed() ? pos_ + rawSize() - 1 - idx : pos_ + idx; }
    /** Helper function that returns the position in the underlying data for a given normalized io-index and symbol type. */
    long indexToPos(long idx = 0, bool outputSymbol = true) const {
		if (nU() == 0) return rawIndexToPos(idx);
		return outputSymbol + ( isReversed() ? 2*((pos_ + rawSize() - 1)/2 - idx) : 2*(pos_/2 + idx) );
	}
    /** Check that the given position `p` in the underlying data is valid, i.e., assert `pos_ <= p < pos_ + rawSize()`, and return `p`. */
    long checkPos(long p) const CHECK(throw(std::out_of_range)) {
        CHECK( if (p < pos_ or p >= pos_ + rawSize()) { throw std::out_of_range("sequence index out of range"); } )
        return p;
    }
    
}; // class Sequence

// The following is some magic to make Sequence objects slice-able and iterable in python
// MARK: python interface stuff
SWIGCODE(%ignore stop_iteration;)
struct stop_iteration {};
#ifdef SWIG
%extend Sequence {
	%typemap(in) PySliceObject* {
		if (!PySlice_Check($input)) { %argument_fail(SWIG_TypeError, "$type", $symname, $argnum); }
		$1 = (PySliceObject *) $input;
	}
	%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) PySliceObject* { $1 = PySlice_Check($input); }
	%feature("python:slot", "nb_nonzero", functype="inquiry") __nonzero__;
	bool __nonzero__() const { return ($self->nO() != 0); }
	bool __bool__() const { return ($self->nO() != 0); }
	%feature("python:slot", "mp_subscript", functype="binaryfunc") __getitem__;
	Sequence __getitem__(PySliceObject *slice) throw (std::invalid_argument CHECK(, std::out_of_range)) {
		Py_ssize_t start, stop, step, sliceLen;
		PySlice_GetIndicesEx(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)$self->rawSize(), &start, &stop, &step, &sliceLen);
		if (step != 1 and step != -1) { throw std::invalid_argument("slice step must be -1 or 1 for tom::Sequence objects."); }
        if (start < 0) { start = tom::NoIndex; }
        if (stop < 0) { stop = tom::NoIndex; }
        return $self->rawSlice(start, stop, step == 1);
    }
    Symbol __getitem__(long i) throw (std::out_of_range) {
        if (i < -$self->rawSize() or i >= $self->rawSize()) { throw std::out_of_range("sequence index out of range"); }
        return $self->rawAt(i);
    }
    %feature("python:slot", "mp_ass_subscript", functype="objobjargproc") __setitem__;
    void __setitem__(long i, Symbol val) throw (std::out_of_range) {
        if ((i < -$self->rawSize()) or (i >= $self->rawSize())) { throw std::out_of_range("sequence index out of range"); }
        $self->rawAt(i, val);
    }
    %feature("python:slot", "tp_iter", functype="getiterfunc") __iter__;
    Sequence __iter__() { return tom::Sequence(*$self); }
    %typemap(throws) tom::stop_iteration {
    (void)$1;
    SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
    SWIG_fail;
    }
	%feature("python:slot", "tp_iternext", functype="iternextfunc") __next__;
    Symbol __next__() throw(tom::stop_iteration) {
        if ($self->rawSize() == 0) { throw tom::stop_iteration(); }
        tom::Symbol ret = $self->rawAt(0); $self->incr_as_python_iterator_only(); return ret;
	}
};
#endif // SWIG

/** A vector of sequences */
typedef std::vector<Sequence> Sequences;

#ifndef SWIG

//IO-functions
std::ostream& Sequence::writeFormattedData(std::ostream& ostream, long maxOut) const {
    if (maxOut == 0) { maxOut = length(); }
    char io_char = isReversed() ? '>' : '<';
    char gap_char = ' ';
    bool theres_more = maxOut < length();
    maxOut = std::min(maxOut, length());
    if (!isIO()) {
        if (maxOut == 0) return ostream;
        for (unsigned long i = 0; i < maxOut - 1; ++i) {
            ostream << rawAt(i) << gap_char;
        }
        ostream << rawAt(maxOut-1);
    } else {
        maxOut *= 2;
        maxOut = std::min(maxOut, rawSize());
        if (maxOut == 0) return ostream;
        bool align = false;
        if ( (isReversed() and !isBackAligned() ) or ( !isReversed() and !isFrontAligned() ) ) {
            ostream << gap_char << io_char;
            align = true;
        }
        for (unsigned long i = 0; i < maxOut-1; ++i) {
            align = not align;
            ostream << rawAt(i) << ( align ? io_char : gap_char );
        }
        ostream << rawAt(maxOut-1);
        if (!align) { ostream << io_char; }
    }
    if (theres_more) { ostream << " ..."; }
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
