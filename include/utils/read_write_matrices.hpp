#ifndef READ_WRITE_MATRICES_HPP
#define READ_WRITE_MATRICES_HPP

// EXTERNAL LIBs ------
#ifdef USE_EIGEN_PKGCONFIG_HEADERS
    #include <eigen3/Eigen/Dense>
    #include <eigen3/Eigen/StdVector>
#else
    #include <Eigen/Dense>
    #include <Eigen/StdVector>
#endif

// STL/std ------------
#include <typeinfo>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>




namespace amp {

template<typename Scalar, typename Derived>
struct rwmatrix
{
    typedef typename Eigen::DenseBase<Derived>::Index iter_t;
    //! appends integer to file name
    static const std::string append_int(const char* filename, const int number);

    //! appends integer and .dat to file name
    static const std::string generic_name(const char* filename, const int number, const bool dat=true);

    static const std::string generic_name(const char* filename, const char* extra, const int number, const bool dat=true);

    //! writes the data in \a mat_in to a file \a filename
    static void write_matrix_to_file(const char* filename, const Eigen::DenseBase<Derived>& mat_in);

    //! writes the data in \a mat_in to a file \a filename
    static void write_matrix_to_file(const char* filename, const Eigen::MapBase<Derived>& mat_in);

    //! writes the data in \a mat_in to a file \a filename
    static void write_matrix_to_file(const char* filename,
                                     const Eigen::Map<Eigen::Matrix<Scalar,-1,-1,Eigen::RowMajor>,0,
                                     Eigen::Stride<-1, -1> >& mat_in);

    //! writes the data in \a mat_in to a file \a filename
    static void write_matrix_to_file(const char* filename,
                                      Eigen::Map<const Derived, 0, Eigen::Stride<-1, -1> >& mat_in);

    //! writes the data in \a mat_in to a file \a filename
    static void write_csl_matrix_to_file(const char* filename, const Eigen::DenseBase<Derived>& mat_in);

    //! writes the data in \a mat_in to a file \a filename
    static void write_csl_matrix_to_file(const char* filename, const Eigen::MapBase<Derived>& mat_in);

    //! writes the data in \a mat_in to a file \a filename
    static void write_csl_matrix_to_file(const char* filename,
                                     const Eigen::Map<Eigen::Matrix<Scalar,-1,-1,Eigen::RowMajor>,0,
                                     Eigen::Stride<-1, -1> >& mat_in);


    //! writes the data in \a mat_in to a file \a filename
    static void write_csl_matrix_to_file(const char* filename,
                                      Eigen::Map<const Derived, 0, Eigen::Stride<-1, -1> >& mat_in);


    //! reads the data from \a filename and writes it to an Eigen matrix \a mat_in. The size
    //! of the data being read in must be known a priori
    static void read_matrix_from_file(const char* filename, Eigen::DenseBase<Derived>& mat_in, const iter_t& r, const iter_t& c);

    static void read_csl_matrix_from_file(const char* filename, Eigen::DenseBase<Derived>& mat_in, const int &r, const int &c);

    //! reads the data from \a filename and writes it to an Eigen matrix \a mat_in. Does not require
    //! a priori information about size, but must read through data 2 times.
    static void size_read_matrix_from_file(const char* filenm, Eigen::DenseBase<Derived>& mat);

    //! reads the data from a one line string and writes it to an Eigen matrix \a mat_in. Does not require
    //! a priori information about size, but must read through data 2 times.
    static void size_read_vector_from_string(const char* string_, Eigen::DenseBase<Derived>& mat);
};



//! append integer and .dat to file name
//!
template<typename Scalar, typename Derived>
const std::string rwmatrix<Scalar, Derived>::
append_int(const char* filename, const int number)
{
    std::stringstream ss_name;
    ss_name << filename << "_" << number;
    return( ss_name.str() );
}




//! append integer and .dat to file name
//!
template<typename Scalar, typename Derived>
const std::string rwmatrix<Scalar, Derived>::
generic_name(const char* filename, const int number, const bool dat)
{
    std::stringstream ss_name;
    ss_name << filename << "_" << number;
    if (dat) {ss_name << ".dat";}
    return( ss_name.str() );
}



//! append integer and .dat to file name
//!
template<typename Scalar, typename Derived>
const std::string rwmatrix<Scalar, Derived>::
generic_name(const char* filename, const char* suffix, const int number, const bool dat)
{
    std::stringstream ss_name;
    ss_name << filename << "_" << suffix << "_" << number;
    if (dat) {ss_name << ".dat";}
    return( ss_name.str() );
}


template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>::
write_matrix_to_file(const char* filename, const Eigen::DenseBase<Derived>& mat_in)
{
    iter_t r = 0;
    iter_t c = 0;
    r = mat_in.rows();
    c = mat_in.cols();

    // setup output file
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::out);

    // print data to file
    for (iter_t i=0; i<r; i++) {
        std::ostringstream line;
        line.setf(std::ios::fixed, std::ios::floatfield);
        line.precision(6);
        for (iter_t j=0; j<c; j++) {
            line << std::setw(12) << mat_in(i,j) << " ";
        }
        std::string strs = line.str();
        output_file << strs << std::endl;
    }
    output_file.close();
}

//! Function to write values to a file
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>
::write_matrix_to_file (const char* filename, const Eigen::MapBase<Derived>& mat_in)
{
    iter_t r = 0;
    iter_t c = 0;
    r = mat_in.rows();
    c = mat_in.cols();

    // setup output file
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::out);

    // print data to file
    for (iter_t i=0; i<r; i++) {
        std::ostringstream line;
        line.setf(std::ios::fixed, std::ios::floatfield);
        line.precision(6);
        for (iter_t j=0; j<c; j++) {
            line << std::setw(12) << mat_in(i,j) << " ";
        }
        std::string strs = line.str();
        output_file << strs << std::endl;
    }
    output_file.close();
}



//! Function to write values to a file
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>
::write_matrix_to_file(const char* filename,
                       Eigen::Map<const Derived, 0, Eigen::Stride<-1, -1> > &mat_in)
{
    iter_t r = 0;
    iter_t c = 0;
    r = mat_in.rows();
    c = mat_in.cols();

    // setup output file
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::out);

    // print data to file
    for (iter_t i=0; i<r; i++) {
        std::ostringstream line;
        line.setf(std::ios::fixed, std::ios::floatfield);
        line.precision(6);
        for (iter_t j=0; j<c; j++) {
            line << std::setw(12) << mat_in(i,j) << " ";
        }
        std::string strs = line.str();
        output_file << strs << std::endl;
    }
    output_file.close();
}




//! Function to write values to a file
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>
:: write_matrix_to_file(const char* filename,
                        const Eigen::Map<Eigen::Matrix<Scalar,-1,-1,Eigen::RowMajor>,0,
                        Eigen::Stride<-1, -1> >& mat_in)
{
    iter_t r = 0;
    iter_t c = 0;
    r = mat_in.rows();
    c = mat_in.cols();

    // setup output file
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::out);

    // print data to file
    for (iter_t i=0; i<r; i++) {
        std::ostringstream line;
        line.setf(std::ios::fixed, std::ios::floatfield);
        line.precision(6);
        for (iter_t j=0; j<c; j++) {
            line << std::setw(12) << mat_in(i,j) << " ";
        }
        std::string strs = line.str();
        output_file << strs << std::endl;
    }
    output_file.close();
}











template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>::
write_csl_matrix_to_file(const char* filename, const Eigen::DenseBase<Derived>& mat_in)
{
    iter_t r = 0;
    iter_t c = 0;
    r = mat_in.rows();
    c = mat_in.cols();

    // setup output file
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::out);

    // print data to file
    for (iter_t i=0; i<r; i++) {
        std::ostringstream line;
        line.setf(std::ios::fixed, std::ios::floatfield);
        line.precision(6);
        for (iter_t j=0; j<c; j++) {
            line << std::setw(12) << mat_in(i,j) << ", ";
        }
        std::string strs = line.str();
        output_file << strs << std::endl;
    }
    output_file.close();
}


//! Function to write values to a file
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>
::write_csl_matrix_to_file (const char* filename, const Eigen::MapBase<Derived>& mat_in)
{
    iter_t r = 0;
    iter_t c = 0;
    r = mat_in.rows();
    c = mat_in.cols();

    // setup output file
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::out);

    // print data to file
    for (iter_t i=0; i<r; i++) {
        std::ostringstream line;
        line.setf(std::ios::fixed, std::ios::floatfield);
        line.precision(6);
        for (iter_t j=0; j<c; j++) {
            line << std::setw(12) << mat_in(i,j) << ", ";
        }
        std::string strs = line.str();
        output_file << strs << std::endl;
    }
    output_file.close();
}



//! Function to write values to a file
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>
::write_csl_matrix_to_file(const char* filename,
                       Eigen::Map<const Derived, 0, Eigen::Stride<-1, -1> > &mat_in)
{
    iter_t r = 0;
    iter_t c = 0;
    r = mat_in.rows();
    c = mat_in.cols();

    // setup output file
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::out);

    // print data to file
    for (iter_t i=0; i<r; i++) {
        std::ostringstream line;
        line.setf(std::ios::fixed, std::ios::floatfield);
        line.precision(6);
        for (iter_t j=0; j<c; j++) {
            line << std::setw(12) << mat_in(i,j) << ", ";
        }
        std::string strs = line.str();
        output_file << strs << std::endl;
    }
    output_file.close();
}




//! Function to write values to a file
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>
:: write_csl_matrix_to_file(const char* filename,
                        const Eigen::Map<Eigen::Matrix<Scalar,-1,-1,Eigen::RowMajor>,0,
                        Eigen::Stride<-1, -1> >& mat_in)
{
    iter_t r = 0;
    iter_t c = 0;
    r = mat_in.rows();
    c = mat_in.cols();

    // setup output file
    std::ofstream output_file;
    output_file.open(filename, std::ofstream::out);

    // print data to file
    for (iter_t i=0; i<r; i++) {
        std::ostringstream line;
        line.setf(std::ios::fixed, std::ios::floatfield);
        line.precision(6);
        for (iter_t j=0; j<c; j++) {
            line << std::setw(12) << mat_in(i,j) << ", ";
        }
        std::string strs = line.str();
        output_file << strs << std::endl;
    }
    output_file.close();
}











//! Read in data from a text file and set values of a pre-allocated matrix object
//! \a r is the number rows, \a c ... columns in \a mat_in
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>::
read_matrix_from_file(const char* filename, Eigen::DenseBase<Derived>& mat_in, const iter_t& r, const iter_t& c)
{
    // setup input file
    std::ifstream input_file(filename);
    Scalar buffer = 0.0;

    // print data to file
    if (input_file.is_open())
    {
        while( !input_file.eof() ) {
            for (iter_t i=0; i<r; i++) {
                std::string line;
                std::istringstream linestream;
                std::getline(input_file, line);
                // check for commented lines
                while(line.compare(0,1,"#")==0) {
                    std::getline(input_file, line);
                }
                if (line.compare(0,4,"#eof")==0)
                    break;
                else if (input_file.eof())
                    break;

                // check for commented lines
                if (line.compare(0,1,"#")==0){
                    continue;
                }
                linestream.str(line);
                for (iter_t j=0; j<c; j++) {
                    std::string elem;
                    std::getline(linestream, elem, ' ');
                    std::stringstream(elem) >> buffer;
                    mat_in(i,j) = buffer;
                }
            }
        }
    }
    input_file.close();
}



//! Read in data from a text file and set values of a pre-allocated matrix object
//! \a r is the number rows, \a c ... columns in \a mat_in
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>
::read_csl_matrix_from_file(const char* filename, Eigen::DenseBase<Derived>& mat_in, const int &r, const int &c)
{
    // setup input file
    Scalar buff = 0.0;
    std::ifstream input_file(filename);

    // print data to file
    if (input_file.is_open())
    {
        for (int i=0; i<r; i++) {
            std::string line;
            for (int j=0; j<c-1; j++) {
                std::getline(input_file, line, ',');
                if( input_file.eof() )
                    break;
                std::stringstream(line) >> buff;
                mat_in(i,j) = buff;
            }
            if( input_file.eof() )
                break;
            std::getline(input_file, line, ',');
            std::stringstream(line) >> buff;
            mat_in(i,c-1) = buff;
        }
    }
    input_file.close();
}


//! Read in data from a text file and set values of a dynamically allocated matrix object
//!
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>
::size_read_matrix_from_file(const char* filenm, Eigen::DenseBase<Derived>& mat)
{
    // setup input file
    iter_t r_ = 0;
    iter_t c_ = 0;
    std::ifstream input_file(filenm);

    // print data to file
    if (input_file.is_open())
    {
        // find out how many rows and columns we're dealing with...
        while( !input_file.eof() ) {
            std::string line;
            std::istringstream linestream;
            std::getline(input_file, line);

            // check for commented lines
            if (line.compare(0,1,"#")==0)
                continue;
            else if (line.compare(0,4,"#eof")==0)
                break;
            else if (input_file.eof())
                break;
            else
                r_++;

            // extract each column from line
            linestream.str(line);
            iter_t ctmp = 0;
            while(!linestream.eof()) {
                std::string elem;
                std::getline(linestream, elem, ' ');
                if (!linestream.eof())
                    ctmp++;
                else if(elem.size()>0)
                    ctmp++;
                else
                    break;
            }
            c_ = ctmp;
        }
    }
    input_file.close();

    // set size of mat and copy values...
    mat = Eigen::DenseBase<Derived>::Zero(r_,c_);
    read_matrix_from_file(filenm,mat,r_,c_);
}




//! Read in data from a text file and set values of a dynamically allocated matrix object
//!
template<typename Scalar, typename Derived>
void rwmatrix<Scalar, Derived>::
size_read_vector_from_string(const char* string_,  Eigen::DenseBase<Derived>& mat)
{
    // setup input file
    iter_t c_ = 0;

    // find out how many columns we're dealing with...
    std::istringstream linestream;
    std::istringstream linestream_;

    // extract each column from line
    linestream.str(string_);
    iter_t ctmp = 0;
    while(!linestream.eof()) {
        std::string elem;
        std::getline(linestream, elem, ' ');
        if (!linestream.eof())
            ctmp++;
        else if(elem.size()>0)
            ctmp++;
        else
            break;
    }
    c_ = ctmp;

    // set size of mat and copy values...
    Scalar buffer = 0.0;
    linestream_.str(string_);
    mat = Eigen::DenseBase<Derived>::Zero(1,c_);
    for (int i=0; i<c_; i++) {
        std::string elem;
        std::getline(linestream_, elem, ' ');
        std::stringstream(elem) >> buffer;
        mat(0,i) = buffer;
    }
}


} //end namespace amp
#endif // READ_WRITE_MATRICES_HPP
