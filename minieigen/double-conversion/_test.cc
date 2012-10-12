
#include"bignum.cc"
#include"bignum-dtoa.cc"
#include"cached-powers.cc"
#include"diy-fp.cc"
#include"double-conversion.cc"
#include"fast-dtoa.cc"
#include"fixed-dtoa.cc"
#include"strtod.cc"

#include<boost/lexical_cast.hpp>

#include<string>

double_conversion::DoubleToStringConverter doubleToString(
	double_conversion::DoubleToStringConverter::NO_FLAGS,
	"inf", /* infinity symbol */
	"nan", /* NaN symbol */
	'e', /*exponent symbol*/
	-5, /* decimal_in_shortest_low: 0.0001, but 0.00001->1e-5 */
	7, /* decimal_in_shortest_high */
	/* the following are irrelevant for the shortest representation */
	6, /* max_leading_padding_zeroes_in_precision_mode */
	6 /* max_trailing_padding_zeroes_in_precision_mode */
);

std::string doubleToShortest(double d){
	char buf[32];
	double_conversion::StringBuilder sb(buf,32);
	doubleToString.ToShortest(d,&sb);
	return std::string(sb.Finalize());
} 

int main(int argc, char** argv){

	for(int i=1; i<argc; i++){
		double d=boost::lexical_cast<double>(argv[i]);
		std::cerr<<doubleToShortest(d)<<std::endl;
	}
}
