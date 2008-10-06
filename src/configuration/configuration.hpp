/*       Created   :  10/03/2008 09:40:20 PM
 *       Last Change: Mon Oct 06 05:00 PM 2008 CEST
 */
#ifndef __CONFIGURATION_HPP__
#define __CONFIGURATION_HPP__

#include <boost/shared_ptr.hpp>

// forward-declare any, options_description
namespace boost{
	class any;
	namespace program_options{
		class options_description;
	}
}

class Configuration{
	// the public interface is doubled by Configuration_Impl,
	// where the real work is done.
	private:
		struct Configuration_Impl;
		boost::shared_ptr<Configuration_Impl> mImpl;
	public:
		Configuration();
		void addModuleOptions(const boost::program_options::options_description& od);
		int parsecfg(int argc, char* argv[]);
		void conflicting_options(const std::string& s, const std::string& t);
		void dependent_options(const std::string& s, const std::string& t);

		boost::any get(const std::string& s);

		// convenience stuff for get
		std::string getString(const std::string& s);
		float getFloat(const std::string& s);
		int   getInt(const std::string& s);
		bool  getBool(const std::string& s);
};

#endif /* __CONFIGURATION_HPP__ */
