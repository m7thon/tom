/**
 * @file   CerealTom.h
 * @author Michael Thon
 *
 * @brief  This file provides some "cerealization" support functions.
 */

#ifndef CEREAL_TOM_H
#define CEREAL_TOM_H

namespace tom {
#define MVAR(AR, T) AR(cereal::make_nvp(#T, T ## _))
#define VAR(AR, T) AR(cereal::make_nvp(#T, T))
#define OMVAR(AR, T) try { MVAR(AR, T); } catch(...) {}
#define OVAR(AR, T) try { VAR(AR, T); } catch(...) {}
#define INSERT_JSON_IO_FUNCTIONS()																			\
	std::string toJSON() const {																					\
		std::stringstream oss;																							\
		{																																		\
			cereal::JSONOutputArchive ar( oss );															\
      save( ar );																												\
		}																																		\
		return oss.str();																										\
	}																																			\
	void fromJSON(const char * string ) {																	\
		{																																		\
			cereal::JSONInputArchive ar( string );														\
			load( ar );																												\
		}																																		\
	}
}

namespace cereal {

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> void
save(Archive & ar, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m) {
	ar.makeArray();
	for (unsigned long i = 0; i < m.rows(); ++i) {
		ar.startNode(); ar.makeArray(); ar.writeName();
		for (unsigned long j = 0; j < m.cols(); ++j) { ar(m(i,j)); }
		ar.finishNode();
	}
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> void
save(Archive & ar, const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m) {
	ar.makeArray();
	for (unsigned long i = 0; i < m.rows(); ++i) {
		ar.startNode(); ar.makeArray(); ar.writeName();
		for (unsigned long j = 0; j < m.cols(); ++j) { ar(m(i,j)); }
		ar.finishNode();
	}
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> void
load(Archive & ar, const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m) {
	size_type rows, cols; ar.loadSize(rows);
	for (unsigned long i = 0; i < rows; ++i) {
		ar.startNode();
		if (i == 0) {
			ar.loadSize(cols);
			const_cast<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>&>(m).resize(rows, cols);
		}
		for (unsigned long j = 0; j < m.cols(); ++j) {
			ar(const_cast<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>&>(m)(i,j));
		}
		ar.finishNode();
	}
}

template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> void
load(Archive & ar, const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m) {
	size_type rows, cols; ar.loadSize(rows);
	for (unsigned long i = 0; i < rows; ++i) {
		ar.startNode();
		if (i == 0) {
			ar.loadSize(cols);
			const_cast<Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>&>(m).resize(rows, cols);
		}
		for (unsigned long j = 0; j < m.cols(); ++j) {
			ar(const_cast<Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>&>(m)(i,j));
		}
		ar.finishNode();
	}
}

} // namespace cereal

#endif // CEREAL_TOM_H
