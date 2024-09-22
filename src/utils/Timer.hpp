#pragma once
#include <chrono>

#ifdef USE_C_TIMER
#include <time.h>
#endif

/*
Class for timers
*/

namespace element_validity {
	#ifdef USE_C_TIMER
	class Timer {
		private:
		volatile double startTime = 0;
		volatile double duration = 0;
		bool running = false;

		public:
		void start() {
			if (!running) {
				startTime = clock();
				running = true;
			}
		}
		void stop() {
			if (running) {
				double stopTime = clock();
				duration = stopTime - startTime;
				running = false;
			}
		}
		void reset() { startTime = 0; duration = 0; running = false; }

		template <typename U>
		double read() const {
			// std::cerr << duration << " " << CLOCKS_PER_SEC << " ";
			const double ds = duration / CLOCKS_PER_SEC;
			// std::cerr << ds << std::endl;
			if constexpr (std::is_same_v<U, std::chrono::seconds>) {
				return ds;
			}
			else if constexpr (std::is_same_v<U, std::chrono::milliseconds>) {
				return ds*1000;
			}
			else if constexpr (std::is_same_v<U, std::chrono::microseconds>) {
				return ds*1000000;
			}
			else if constexpr (std::is_same_v<U, std::chrono::nanoseconds>) {
				return ds*1000000000;
			}
			else {
				throw std::runtime_error("Unknown time unit");
			}
		}
	};
	#else
	class Timer {
		private:
		using myClock = std::chrono::steady_clock;
		myClock::time_point startTime;
		myClock::duration duration = myClock::duration::zero();
		bool running = false;

		public:
		void start() {
			if (!running) {
				startTime = myClock::now();
				running = true;
			}
		}
		void stop() {
			if (running) {
				duration = myClock::now() - startTime;
				running = false;
			}
		}
		void reset() { duration = myClock::duration::zero(); running = false; }


		template <typename U>
		typename U::rep read() const {
			return std::chrono::duration_cast<U>(duration).count();
		}
	};
	#endif
}