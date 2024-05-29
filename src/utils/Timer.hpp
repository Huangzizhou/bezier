#pragma once
#include <chrono>

/*
Class for timers
*/

namespace element_validity {
	using chrono_t = double;

	class Timer {
		private:
		using clock = std::chrono::high_resolution_clock;
		clock::time_point startTime;
		clock::duration duration = clock::duration::zero();
		bool running = false;

		public:
		void start() { startTime = clock::now(); running = true; }
		void stop() { duration += clock::now() - startTime; running = false; }
		void reset() { duration = clock::duration::zero(); }

		// template <TimeUnitType U>
		template <typename U>
		typename U::rep read() const {
			const auto dur =
				running ? duration + clock::now() - startTime : duration;
			return std::chrono::duration_cast<U>(dur).count();
		}
	};
}