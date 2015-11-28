namespace tom {

SWIGCODE(%feature("director") StopCondition;)

/**
 * This class allows to control the stopping condition for iterative algorithms, and at the same time provides a callback mechanism between iterations.
 */
class StopCondition {
public:
    int maxIterations_; /**< the maximum number of iterations to perform */
    double relativeImprovementThreshold_; /**< stop the iteration if the relative improvement falls below this threshold */
    double absoluteImprovementThreshold_; /**< stop the iteration if the absolute improvement falls below this threshold */
    int iteration_ = 0; /**< the current iteration */
    double lastValue_ = -std::numeric_limits<double>::quiet_NaN(); /**< the most recent "value" achieved by the iteration */

    /** Construct a `StopCondition` while setting the key parameters.
     */
    StopCondition(int maxIterations = 100, double relativeImprovementThreshold = 1e-7, double absoluteImprovementThreshold = 1e-12) :
            maxIterations_(maxIterations), relativeImprovementThreshold_(relativeImprovementThreshold),
            absoluteImprovementThreshold_(absoluteImprovementThreshold) { }

    /** Return `true` if the iteration should be terminated. This method gets called before starting a new iteration. The termination condition can be customized by inheriting from `StopCondition` and overwriting this method (even from Python).
     *
     * By default, the iteration is terminated if either of the following conditions is met:
     *
     * - the `iteration_` count reaches `maxIterations_`
     * - the relative improvement |(`currentValue` - `lastValue_`) / `lastValue_`| is less than the `relativeImprovementThreshold_`
     * - the absolute improvement |(`currentValue` - `lastValue_`)| is less than the `absoluteImprovementThreshold_`
     */
    virtual bool stop(double currentValue) {
        return (std::fabs((currentValue - lastValue_) / lastValue_) < relativeImprovementThreshold_ or
                std::fabs(currentValue - lastValue_) < absoluteImprovementThreshold_ or
                iteration_ >= maxIterations_);
    }

    /** This method gets called before every iteration and also right before termination. It is a no-op by default, but can be used as a callback between iterations by inheriting from `StopCondition` and overwriting this method (even from Python).
     */
    virtual void callback() { }

    /** This method gets called with the new "value" `currentValue` before every iteration and serves several purposes:
     *
     * - it calls `.stop(currentValue)` to decide if the iteration should terminate
     * - the `.lastValue_` is updated to the `currentValue`
     * - it calls `.callback()` to allow for callbacks
     * - it updates the `iteration_` count
     *
     * Note that `inf` or `nan` may be passed before the first iteration if no "value" is available yet. */
    bool operator()(double currentValue = std::numeric_limits<double>::quiet_NaN()) {
        bool stop = this->stop(currentValue);
        lastValue_ = currentValue;
        this->callback();
        if (not stop) { iteration_++; }
        return stop;
    }

    /** Reset the `iteration_` count and the `lastValue_`. */
    void reset() {
        iteration_ = 0;
        lastValue_ = std::numeric_limits<double>::quiet_NaN();
    }

    virtual ~StopCondition() { };
};

} // namespace tom