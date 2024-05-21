# Reproducing results by Cochrane & Piazzesi (2005)

The research of Cochrane & Piazzesi (2005) focuses on examining the changing patterns of expected excess bond returns over time. To investigate this, they perform regressions where we analyze the relationship between one-year excess returns and initial forward rates. The findings reveal that a single factor, represented by a tent-shaped linear combination of forward rates, effectively predicts excess returns across one-year to five-year maturity bonds.

To replicate the analysis of Cochrane & Piazzesi (2005), we begin by calculating the log-yield, denoted as $y_t^{(n)}$, which is defined as $-\frac{1}{n}p_t^{(n)}$, where $p_t^{(n)}$ represents the log price of an n-year discount bond at time t. We then compute the log forward rates, denoted as $f_t^{(n)}$, given by $p_t^{(n-1)} - p_t^{(n)}$. Additionally, we need to obtain the log holding period, denoted as $r_{t+1}^{(n)}$, defined as $p_{t+1}^{(n-1)} - p_t^{(n)}$ where $p_{t+1}^{(n-1)}$ is the log price lagged by 12 months. With these calculations in place, we can determine the excess log returns, denoted as $rx_{t+1}^{(n)}$, which are obtained as $r_{t+1}^{(n)} - y_t^{(1)}$, where $y_t^{(1)}$ represents the log-yield in the first period.

Next, we can compute the excess return forecasts using the unrestricted model expressed as follows:

$$
rx_{t+1}^{(n)} = \beta_0^{(n)} + \beta_1^{(n)}y_t^{(n)} + \beta_2 f_t^{(2)} + \beta_3 f_t^{(3)} + \beta_4 f_t^{(4)} + \beta_5 f_t^{(5)} + \epsilon_{t+1}^{(n)}
$$

Results of this regression are displayed in Figure 1.

To impose restrictions on the model, we need to derive the $\gamma$ values by fitting the following model of the average excess return on all forward rates:

$$
\frac{1}{4}\sum_{n=2}^5 rx_{t+1}^{(n)} = \gamma_0 + \gamma_1 y_t^{(n)} + \gamma_2 f_t^{(2)} + \gamma_3 f_t^{(3)} + \gamma_4 f_t^{(4)} + \gamma_5 f_t^{(5)} + \bar \epsilon_{t+1} 
$$

$$
\bar{rx}{t+1} = \gamma^T f_t + \bar \epsilon{t+1}
$$

Results of this regression are displayed in Panel A of Table 1.

Subsequently, we can estimate the single factor restricted model by conducting four regressions for $n={2,3,4,5}$:

$$
rx_{t+1}^{(n)} = b_n (\gamma^T f_t) + \epsilon_{t+1}^{(n)}
$$

Results of this regression are displayed in Panel B of Table 1 and Figure 1.

Lastly, the final model reproduced in the analysis is the Fama & Bliss (1987) excess return model, given by:

$$
rx_{t+1}^{(n)} = \alpha + \beta(f_t^{(n)}- y_t^{(1)}) + \epsilon_{t+1}^{(n)}
$$

Results of this regression are displayed in Table 2.

Looking at the results of the replicated analysis, we see similar results with Fama-Bliss data from January 1964 to December 2003 as found by Cochrane & Piazzesi (2005). Figure 1 suggests a clear tent-shaped coefficients of one-year excess returns on forward rated, Model (1) and Model (4). This holds for both, unrestricted and restricted models.

<figcaption align = "center"><b>Table 1: Estimates of the single-factor model using the Fama-Bliss yields.</b></figcaption> 
<br />

<table class="table" style="font-size: 10px; ">
 <thead>
<tr><th style="border-bottom:hidden;padding-bottom:0; padding-left:5px;padding-right:5px;text-align: left; font-weight: bold; " colspan="9"><div style>A. Estimates of the return-forecasting factor</div></th></tr>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> $\gamma_0$ </th>
   <th style="text-align:right;"> $\gamma_1$ </th>
   <th style="text-align:right;"> $\gamma_2$ </th>
   <th style="text-align:right;"> $\gamma_3$ </th>
   <th style="text-align:right;"> $\gamma_4$ </th>
   <th style="text-align:right;"> $\gamma_5$ </th>
   <th style="text-align:right;"> $R^2$ </th>
   <th style="text-align:right;"> $\chi^2(5)$ </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;width: 4cm; "> OLS estimates </td>
   <td style="text-align:right;"> -1.19 </td>
   <td style="text-align:right;"> -1.29 </td>
   <td style="text-align:right;"> -0.34 </td>
   <td style="text-align:right;"> 1.73 </td>
   <td style="text-align:right;"> 1.19 </td>
   <td style="text-align:right;"> -1.11 </td>
   <td style="text-align:right;">  </td>
   <td style="text-align:right;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 4cm; "> NW, 18 lags </td>
   <td style="text-align:right;"> 1.88 </td>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:right;"> 0.82 </td>
   <td style="text-align:right;"> 0.70 </td>
   <td style="text-align:right;"> 0.51 </td>
   <td style="text-align:right;"> 0.56 </td>
   <td style="text-align:right;"> 0.2 </td>
   <td style="text-align:right;"> 22.51 </td>
  </tr>
</tbody>
</table>

In Panel A of Table 1, the $\gamma$ values for Model (2) closely align with the calculations made by Cochrane & Piazzesi (2005), exhibiting a similar pattern. However, our analysis reveals larger Newey-West standard errors compared to the original paper. Additionally, the $\chi^2(5)$ statistic is smaller in our analysis. The $R^2$ value of 34% remains consistent with the findings of Cochrane & Piazzesi (2005).

Moving to Panel B of Table 2, the coefficients and significant Large T values obtained for the unrestricted Model (1) are nearly identical to those proposed by Cochrane & Piazzesi (2005). The $R^2$ values for the unrestricted Model (1) are slightly higher than reported in the original paper, while being slightly lower for the restricted Model (4).

Table 2 suggests that our results are again in line findings of Cochrane & Piazzesi (2005). We see that the model fit, based on $R^2$ values, is worse in the Fama & Bliss (1987) model than in the single-factor model propose by Cochrane & Piazzesi (2005).

# Check robustness

## Reproduce Cochrane & Piazzesi (2005) using Gürkaynak et al. (2007) yields

To verify the robustness of the results of Cochrane & Piazzesi (2005) we use Gürkaynak-Sack-Wright (GSW) yields. The GSW yields, while having their own limitations present more information than the Fama-Bliss yields. GSW provide daily yield estimates and estimates of the longest available maturities. Also, the data is being updated on a quarterly basis. 

GSW use a parametric method to estimate smooth yield curves. For their estimation GSW exclude among others securities with less than three months to maturity and all Treasury bills, which leads to pricing mismatches on the short end of the maturity spectrum - a deficiency that GSW acknowledge in their paper.

However, in our analysis we only utilize the same maturities as already provided by Fama & Bliss (1987).

When looking at Figure 2, one can immediately observe that the tent shape of the coefficients is lost. Nonetheless, a certain pattern in the coefficients can still be identified. The pattern is more pronounced in the restricted (1) model than in the unrestricted model (4).

Similarly to the analysis with the FB data, also in Table 3 Panel A the $\bf{\gamma}$ from Model (2) reflects the general pattern found in Figure 2. Using GSW data however yields substantially higher Newey-West corrected standard errors and also a $\chi^2(5)$ statistic of 10.46 which is below the critical value at significance level 0.05 for a one sided test (11.07). The 15\% $R^2$ is also well below the 34\% found using the FB data.

In Panel B we can see that the $b_n$ increases again with maturity and the standard errors are still quite low. But unlike with the FB data the $b_n$ are generally lower and the $R^2$ is not as high either. Also, the $R^2$ are higher in the unrestricted model than in the restricted one. 

Looking at Table 4 we get similar results to Fama & Bliss (1987) also with the GSW data. Yet again, the $R^2$ values are worse in the Fama & Bliss (1987) model than in the single-factor model proposed by Cochrane & Piazzesi (2005).

Overall, one can say that the model by Cochrane & Piazzesi (2005) improves the FB model also with the GSW data - albeit less drastically than with the FB data. That said, the expectation hypothesis can still be rejected.


## Reproduce Cochrane & Piazzesi (2005) using recent Fama & Bliss (1987) yields 

We again reproduce the analysis proposed by Cochrane & Piazzesi (2005) using more recent Fama-Bliss yield data from January 1964 to December 2022.

Based on Figure 3, it appears that the conclusions drawn by Cochrane & Piazzesi (2005) may not hold for more recent Famma-Bliss data. The observed coefficients of one-year excess returns on forward rates do not exhibit a distinct tent shape, particularly for forward rates spanning years three to five. Additionally, the coefficients appear to be smaller in magnitude when compared to those depicted in Figure 1.

In Panel A of Table 5, it is evident that the absolute values of the $\gamma$ coefficients are smaller compared to those obtained by Cochrane & Piazzesi (2005). However, the Newey-West standard errors demonstrate greater similarity to the original paper, as indicated by a significantly lower $\chi^2$ statistic of 22.51. The $R^2$ value of 20% represents a decrease of 15% compared to the original paper.

Moving on to Panel B of Table 5, the coefficients and corresponding Large T values appear to align more closely with the findings of Cochrane & Piazzesi (2005). However, when utilizing a larger dataset, the $R^2$ values of both the restricted and unrestricted models experience a decrease of approximately 10%-12%.


According to Table 6, the Famma-Bliss regression conducted on more recent data exhibits lower beta coefficients and lower $R^2$ values. Additionally, the Newey-West standard errors appear to be relatively consistent with those obtained in Table 2.


## Perform an out-of-sample validation using recent Fama & Bliss (1987) yields 

In order to evaluate out-of-sample performs of all the different models used in this analysis, we need to fit all the models once again using the data from January 1964 to $t-1$ where $t_0$ is January 2004. This means that in the first iteration, the models trained based on data from  January 1964 to December 2003 and these models are used to obtain predictions for January 2004. With this, we keep increasing t by 1 month until we can predict the most recent observation which is December 2022. 

To evaluate out-of-sample performance, we calculate the $R^2$ for out-of sample performance suggested by Campbell & Thompson (2007).

$$
R_{OS}^2 = 1 - {\sum_{t=1}^T (r_t - \hat r_t)^2}/{\sum_{t=1}^T (r_t - \bar r_t)^2}
$$

where $r_t$ is the actual observation of excess return, $\hat r_t$ is the predicted excess returns and $\bar r_t$ is the historical average excess return calculated up until $t-1$.

$R_{OS}^2$ allows us to compare out-of-sample $R^2$ with the in-sample $R^2$ statistic. Positive $R_{OS}^2$ value means that "the predictive regression has lower average mean-squared prediction error than the historical average return" (Campbell & Thompson, 2007), p.1515) whereas the opposite is true from negative $R_{OS}^2$.

Table 7 reveals that all models proposed by Cochrane & Piazzesi (2005) exhibit negative $R_{OS}^2$ values, indicating that the historical average excess return yields smaller predictive errors compared to the fitted models. On the other hand, the Fama & Bliss (1987) model demonstrates positive $R_{OS}^2$ values suggesting that the model performs slightly better on out-of-sample data compared to the historical average.

# Discussion

## Main messages of the papers

The government bond market is an important asset class as it has the second-highest market capitalization among the major asset classes (Randl & Weiss, 2023, p. 1). Similar to the assignment on portfolio sorts in the equity market context, we are interested in factors which can predict future bond risk premia. Our analysis will be based on the relationship between short- and long-term bond yields. The fundamental price relation covering bond yields was majorly influenced by Lutz (1940) who argue that the anticipated gain obtained by maintaining a long-term investment in the money or capital market until it reaches maturity is equivalent to the expected return achieved by continuously renewing a sequence of short-term investments that collectively have the same overall maturity as the long-term investment. This implies the following three properties:

a) Long maturity yield = average of expected future short rates 
b) Forward rate = expected future spot rate
c) Expected holding period returns should be equal across maturities 

In practice this means that, the expected 1-year return on a bond is the same whether you invest in a bond with maturity of one year or invest in a bond with maturity of n-years and sell it after the first year. By taking the theory strictly it would mean that the yield, therefore, the price of bonds is fully explained by the current short term yield and the future expected yields. In this case returns are driven by random walks, and thus, cannot be predicted.Given an efficient market with rational investors and the absence of arbitrage opportunities the expected price in the future should equal the current forward rate. However, in the following sections we will explain that this does not hold in practice, one can expect higher returns for holding long maturity bonds instead of rolling short term bonds.


**Fama & Bliss (1987)**  

With the regression of each excess return against the corresponding forward spread, Fama & Bliss (1987) where among the first to reject the expectation hypothesis. They observed treasury bonds and their monthly prices of zero yield bonds of maturity one to five years, from 1964 to 1985. The forward spread is denoted as the difference between the forward rate of maturity n and the yield on the 1-year bonds. 

Following the expectation hypothesis, where $\alpha$ can be different than zero, $\beta$ must be equal to zero, to allow for constant expected excess returns. However, the results of Fama & Bliss (1987) show that $\beta$ is not equal to zero (see Table 2) and hence the expectation hypothesis can not hold. This analysis has become the "classic evidence against the expectations hypothesis in long-term bonds" (Cochrane & Piazzesi, 2005, p.145) since the regressions results have held up until today since their publication, unlike many other anomalies. 


**Cochrane & Piazzesi (2005)** 

As in the previous analysis, Cochrane & Piazzesi (2005) will use monthly treasury bond prices but this time for a greater period with data ranging between 1964 and 2003. The major difference between the previous approach and Cochrane & Piazzesi (2005) is that while before the estimates were based on only one forward spread, we now consider the spread of five forward rates. 

The analysis starts with the unrestricted model which simply regresses five forwards rates on excess returns on the n-year bond. A tent shape is revealed if we are graphing the slope coefficients, omitting the constants, as a function of the maturity n. Under the null that all estimated coefficients are jointly greater than zero, a rejection of the null implies that the coefficients are different from zero. This again implies that the expectation hypothesis is rejected, as the expectations hypothesis only will accept constant or zero excess returns, i.e. the value of the constant coefficient.

In the following, Cochrane & Piazzesi (2005) progress even further to set up a restricted model using only a single factor which predicts returns on bonds of any maturity. This model is estimated in two steps. First $\gamma$ is estimated by regressing the forward vector on the average excess returns. Second the $b_n$ is estimated for each of the four excess returns in the dataset using Model 5.

As we can see in Figure 1 the tent shape reveals even under the restricted single-factor model. Hence, the main message of the analysis is that with the single factor model we can not only reject the expectation hypothesis but also predict returns on bonds of different maturities with a single factor. However, Cochrane & Piazzesi (2005) need to reject the single factor model even as it produces close parameters to the original model (p.156). This is provided by a General Method of Movements specification which not only says that jointly, the unrestricted and the restricted models are not equal but also rejects the single-factor model (Cochrane & Piazzesi, 2005, p.156-157).

Additionally, the tenet shape distorts with recent data. The authors attribute this to business cycles. However, the similarity across maturities stay more or less intact, meaning that the general idea of the model still works. This claim is reinforced by the inclusion of recent data although with much lower R-squares, with around 10 percentage point drop. Moreover, since the pattern is not as orderly, the discrepancy between the 1 factor model and the inclusion of all 5 forward rates separately is getting wider. Since 2008 came with a big structural change with very low yields, a rethinking of the model might be apt with for instance the inclusion of longer maturities.

**Possible improvements:**

The restricted models $R^2$ is only 0.35, however the unrestricted model has an $R^2$ of 44 percent.
The model might be limited to US-Treasury data only. Testing for different countries/currencies needed to verify. 
To generally improve the return predictability, one could add those features to the model which often offsets the correct predictions. In practice bond returns are sensitive to monetary policy, inflation prospects and model uncertainty. To respect these sensitivities we introduce an incomplete list:

- Introduction of time variation in the regression parameters.
-	Introduction of economic variables which have forecasting power like inflation or employment.
-	Combining the previous suggestions with forecast combination methods. Starting with simple equal-weighted averages of different model specifications more factors can be considered.


## Interpret the results or your analysis

From an investment point of view, by the relatively high explanatory power, the model seems consistent enough to identify underpriced bond yields and make gains on the long term. Simple strategies based on the model suggest:

- When the yield spread between long-term and short-term bonds widens, it implies that the market foresees an upcoming rise in short-term interest rates. In such a scenario, allocating a larger portion of your portfolio towards short-term bonds could lead to higher returns.
- On the other hand, when the yield spread narrows, it signals that the market expects a decrease in short-term interest rates. In this situation, one may want to consider expanding their investment in long-term bonds to take advantage of potentially higher yields in comparison to future lower short-term interest rates.
These are only very simplistic examples, however, what is important to point out is that in the CP model introduces a risk factor, which drives future excess returns. Therefore, one achieves higher returns by exposing their portfolio to more risk (no free-lunch), which in case is the variability of interest rates and other macroeconomic variables driving the forward rates (Brooks & Moskowitz, 2017). Thornton & Valente (2012) find that forward spread predictors, when used to guide the investment decisions of an investor with mean-variance preferences, do not lead to higher out-of-sample Sharpe ratios.

Although empirical studies have demonstrated statistical evidence supporting the predictability of bond returns, there is currently limited proof to suggest that this predictability could have been effectively utilized in real-time to enhance investors' economic utility. Similar to what we described, with the inclusion of recent data, Kessler & Scherer (2009) with an international sample also suggest that the unrestricted model has a reasonable forecasting power for future bond returns. On the other hand, the restricted model does not perform as well. They find that accounting for time varying parameters can lead to more accurate forecasts.

**Bond Issuers:**

- *Timing of Issuance:* The CP model's predictions of future changes in short-term interest rates can be valuable for bond issuers in timing their debt issuances. If the model suggests that short-term interest rates are expected to increase, issuers may consider issuing bonds with longer maturities to lock in lower borrowing costs before rates rise.  
- *Yield Curve Management:* The model's insights into the yield spread between long-term and short-term bonds can help issuers manage the shape of the yield curve. By understanding how changes in short-term interest rates affect the yield spread, issuers can make informed decisions on issuing different maturities to achieve their desired yield curve shape.

**Regulators:**

- *Monetary Policy Decisions:* The CP model's ability to predict changes in short-term interest rates can be useful for central banks and regulators in making monetary policy decisions. If the model indicates that the yield spread is widening and short-term interest rates are expected to rise, central banks may consider tightening monetary policy to control inflationary pressures. Conversely, a narrowing yield spread might suggest the need for looser monetary policy to stimulate economic growth.  
- *Financial Stability:* The CP model's insights can help regulators assess financial stability and systemic risks. By monitoring the yield spread and its relationship with short-term interest rates, regulators can gain insights into the market's expectations and potential vulnerabilities. A widening yield spread, for example, could signal market concerns about future economic conditions or liquidity risks.

# References
