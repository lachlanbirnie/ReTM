# Relative Transfer Matrix (ReTM)
Relative Transfer Matrix (a multi-source Relative Transfer Function).

Work in progress ...

1. [Theory of Relative Transfer Function (ReTF)](#theory-of-relative-transfer-function)
    1. [Ideal ReTF](#Ideal-ReTF)
    1. [Obtaining ReTF in Practice](#obtaining-retf-in-practice)
2. [Theory of Relative Transfer Matrix (ReTM)](#theory-of-relative-transfer-matrix)
    1. [Ideal ReTM](#Ideal-ReTM)
    1. [Obtaining ReTM in Practice](#obtaining-retm-in-practice)
3. [Application: Denoising with ReTF/ReTF](#application-denoising-with-retf--retm)

## Theory of Relative Transfer Function
Consider a single sound source with signal $s(t)$. Two microphones $\{1,2\}$ whose recorded signals are $p_{1}(t)$, $p_{2}(t)$. The source and microphones are in a room such that the recorded signals are decribed by the room's impulse responses (RIRs) $h_{\{1,2\}}^{s}(t)$, 
$$
p_{1}(t) = h_{1}^{s}(t) * s(t), \qquad p_{2}(t) = h_{2}^{s}(t) * s(t).
$$

In the short-time frequency domain the recorded signals can be desribed as*
$$
P_{1}(\tau,f) = H_{1}^{s}(f) \times S(\tau,f), \qquad P_{2}(\tau,f) = H_{2}^{s}(f) \times S(\tau,f),
$$
where:
- $(\tau, f)$ denote (time, frequency) bin index.
- $P$, $H$, and $S$ are frequency domain equivalents of $p$, $h$, and $s$.
- \* (so long as the STFT window length is greater than the RIR length, and the NFFT size of the STFT operation is zero-padded to avoid circular aliasing such that NFFT > RIR_length + window_length -1).

### Ideal ReTF
Ideally, the ReTF describes the coupling of two microphones in response to the single sound source. In a stationary environment the ReTF is given by
$$
\begin{equation}
R_{12}^{(s)}(\tau,f) = \frac{H_{1}^{s}(\tau,f)}{H_{2}^{s}(\tau,f)}.
\end{equation}
$$
We note that the ReTF is source-dependent, denoted by the $\cdot^{(s)}$ superscript, and that each source-microphone pairing would have its own unique ReTF.

The ReTF can describe what the recording of microphone 1 will be given the observation of microphone 2's recording of the sound source, described as
$$
\begin{split}
P_{1}(\tau,f) & \equiv R_{12}^{(s)}(\tau,f) \times P_{2}(\tau,f) \\
& \equiv \frac{H_{1}^{s}(\tau,f)}{H_{2}^{s}(\tau,f)} \times H_{2}^{s}(f) \times S(\tau,f) \\
& \equiv H_{1}^{s}(f) \times S(\tau,f).
\end{split}
$$

### Obtaining ReTF in Practice
In practice, the (impulse respones) transfer functions $H_{\{1,2\}}^{s}(\tau, f)$ are unknown for us to get the ReTF with (1). However, the unique property of the ReTF is that it can be estimated from the recorded signals $P_{1}(\tau,f)$ and $P_{2}(\tau,f)$ directly. This sub-section details how:

#### Cross-Power Spectral Density
The cross power spectral density between the microphone signals is given by
$$
\Rho_{1,2}(\tau,f) = P_{1}(\tau,f) \times P^{*}_{2}(\tau,f)
$$
I also refer to the auto-correlation PSD as cross power spectral density of a microphone with itself
$$
\Rho_{2,2}(\tau,f) = P_{2}(\tau,f) \times P^{*}_{2}(\tau,f)
$$

#### Spectral Covariance
Time averaging of the CPSDs to get the spectral covariance
$$
\mathcal{P}_{1,2}(f) = \frac{1}{T} \sum_{t=1}^{T} \Rho_{1,2}(\tau-t,f), 
\qquad 
\mathcal{P}_{2,2}(f) = \frac{1}{T} \sum_{t=1}^{T} \Rho_{2,2}(\tau-t,f).
$$

#### Estimated Relative Transfer Function
The ReTF is estimated by the ratio of spectral covariances
$$
R_{12}^{(s)}(f) \approx \frac{\mathcal{P}_{1,2}(f)}{\mathcal{P}_{2,2}(f)}
$$
Again, such that microphone-1 can be determined from microphone-2
$$
P_{1}(\tau,f) \approx R_{12}^{(s)}(f) \times P_{2}(\tau,f)
$$

## Theory of Relative Transfer Matrix
Now consider multiple sound sources indexed by $\ell = [1, ..., L]$ with signals $s_{\ell}(t)$. In the STFT-domain the set of source signals is given by
$$
\textbf{S}(\tau,f) = \begin{bmatrix} S_{1}(\tau,f) \\ \vdots \\ S_{L}(\tau,f) \end{bmatrix}
$$

Consider two multi-channel microphones, or two groups of microphones:
- The first group, lets call the target group, is denoted "group $A$" with a total of $A$ microphones and recorded signals
$$
\textbf{A}(\tau,f) = \begin{bmatrix} P_{1}(\tau,f) \\ \vdots \\ P_{A}(\tau,f) \end{bmatrix}
$$
- The second group, lets call the reference group, is denoted "group $B$" with a total of $B$ microphones and recorded signals
$$
\textbf{B}(\tau,f) = \begin{bmatrix} P_{1}(\tau,f) \\ \vdots \\ P_{B}(\tau,f) \end{bmatrix}
$$
- Note that $A$ does not need to equal $B$, the two groups can have a different number of microphones.

In a reverberant environment described by transfer functions between the sources and microphones, the received signals for each group of microphones is given by
$$
\textbf{A}(\tau,f) = \textbf{H}_{A}(f) \textbf{S}(\tau,f), \qquad 
\textbf{B}(\tau,f) = \textbf{H}_{B}(f) \textbf{S}(\tau,f)
$$
where
$$
\textbf{H}_{A}(f) = 
\begin{bmatrix} 
H_{1}^{1}(f) & \cdots & H_{1}^{L}(f) \\
\vdots & \ddots & \vdots \\
H_{A}^{1}(f) & \cdots & H_{A}^{L}(f) \\
\end{bmatrix}
$$
and $H_{A}^{L}(\tau,f)$ denotes the room transfer function between the $L$-th source and $A$-th micropohone. $\textbf{H}_{B}(\tau,f)$ is defined similarly.

### Ideal ReTM
Ideally, the ReTM describes the same spatial coupling as the ReTF. That is, the ReTM describes the coupling between two groups of microphones in reponse to a set of sound sources. In a stationary environment the ReTM is given by 
$$
\begin{equation}
\textbf{R}_{AB}(f) = \textbf{H}_{A}(f) \textbf{H}_{B}^{-1}(f).
\end{equation}
$$
In this way, the ReTM can describe the recordings of microphones in group A from an observation of the recordings of group B. 
$$
\begin{split}
\textbf{A}(\tau,f) & = \textbf{R}_{AB}(f) \textbf{B}(\tau,f) \\
& = \textbf{H}_{A}(f) \textbf{H}_{B}^{-1}(f) \textbf{H}_{B}(f) \textbf{S}(\tau,f) \\
& = \textbf{H}_{A}(f) \textbf{S}(\tau,f)
\end{split}
$$

### Obtaining ReTM in Practice
In practice, the room's impulse responses are unknown, and thus the matrices $\textbf{H}_{\{A,B\}}(f)$ are also unknown for getting the ideal ReTM in (2). However, like the ReTF, the ReTM can also be estimated from the recorded signals directly.

#### Cross-Power Spectral Density
Cross-power spectral density between microphone groups $A$ and $B$ is given by
$$
\Rho_{B,A}(\tau,f) = \textbf{B}(\tau,f) \textbf{A}^{*}(\tau,f)
$$
where, $[\cdot]^{*}$ denotes conjugate transpose.

And the cross-power spectral density between microphone within group $A$ is given by
$$
\Rho_{A,A}(\tau,f) = \textbf{A}(\tau,f) \textbf{A}^{*}(\tau,f)
$$

#### Spatial Covariance
Time average CPSDs to get spatial coveraiances
$$
\mathcal{P}_{B,A}(f) = \frac{1}{T} \sum_{t=1}^{T} \Rho_{B,A}(\tau-t,f)
$$
$$
\mathcal{P}_{A,A}(f) = \frac{1}{T} \sum_{t=1}^{T} \Rho_{A,A}(\tau-t,f)
$$

#### Pseudoinverse
Find the pseudoinverse of $\mathcal{P}_{B,A}(f)$, for example
$$
\mathcal{P}_{B,A}^{\dagger}(f) = \left( \mathcal{P}_{B,A}^{*}(f) \mathcal{P}_{B,A}(f) \right)^{-1} \mathcal{P}_{B,A}^{*}(f)
$$

#### Estimate Relative Transfer Matrix
The ReTM can then be estimated by
$$
\textbf{R}_{AB}(\tau,f) \approx \mathcal{P}_{A,A}(f) \times \mathcal{P}_{B,A}^{\dagger}(f)
$$
And can be used for 
$$
\textbf{A}(\tau,f) \approx \textbf{R}_{AB}(f) \textbf{B}(\tau,f)
$$
which can be thought of estimating the target recordings (in group $A$) from the reference recordings (in group $B$).

## Application: Denoising with ReTF / ReTM
This section gives an example of using the ReTF/ReTM to denoise a noisey-speech signal. The following will be formulated for the use of ReTF but it can be easily extended to use ReTM.

For this to work, some assumptions are needed:
- We assume that the microphones do not move during the recording.
- We assume that there is only one noise source for the ReTF case (or a few sources for ReTM case),
- That the noise's signal is mostly continuously active over time,
- That the noise source does not move during the recording.
- We assume that there is one (or more) speech source,
- That the speech is mostly inactive over time, 
- Such that on average over time the recording is mostly noise, this way the ReTF can be estimated from the average of the recording as the averaging will remove most of the speech.

### Problem Formulation
Consider:
- Two microphones, denoted $\{1,2\}$, with recorded signals $p_{1}(t)$, $p_{2}(t)$.
- A noise source, denoted by "$n$" with continuous signal $n(t)$.
- A speech source, denoted by "$s$" with a mostly inactive signal $s(t)$. 
- A reverberant environment described by impulse response between the sources and microphones, $h_{1}^{s}(t)$, $h_{2}^{n}(t)$, ... etc. 

In the STFT-domain, the microphone's recorded siganls are given by
$$
\begin{split}
P_{1}(\tau,f) & = H_{1}^{n}(f) N(\tau,t) + H_{1}^{s}(f) S(\tau,f) \\
P_{2}(\tau,f) & = H_{2}^{n}(f) N(\tau,t) + H_{2}^{s}(f) S(\tau,f)
\end{split}
$$
The goal of ReTF denosing is to enhance an estimate of $S(\tau,f)$ by removing the noise $N(\tau,f)$ from the recordings.

### ReTF Denoising
In the following, we formulate microphone-$1$ as the *target microphone*, and microphone-$2$ as the *reference microphone*. Our goal is to remove the noise from the target microphone's recording by estimating the noise from the reference microphone.

Consider first that we know the ReTF of the noise source $R_{12}^{(n)}(f)$, which idealy would be given by
$$
R_{12}^{(n)}(f) = \frac{H_{1}^{n}(f)}{H_{2}^{n}(f)}
$$

Filtering the reference microphone ($2$) with the ReTF:
$$
\begin{split}
R_{12}^{(n)}(f) P_{2}(\tau,f) & = R_{12}^{(n)}(f) H_{2}^{n}(f) N(\tau,t) + R_{12}^{(n)}(f) H_{2}^{s}(f) S(\tau,f) \\
& = \underbrace{H_{1}^{n}(f) N(\tau,t)}_{\text{estimate of noise}} + R_{12}^{(n)}(f) H_{2}^{s}(f) S(\tau,f)
\end{split}
$$
gives us an estimate of the noise component at the target microphone.

Subtracting this estimate from the target microphone's recording $P_{1}(\tau,f)$ gives:
$$
\begin{split}
P_{1}(\tau,f) - R_{12}^{(n)}(f) P_{2}(\tau,f) = \ & H_{1}^{n}(f) N(\tau,t) + H_{1}^{s}(f) S(\tau,f) 
\\ & - H_{1}^{n}(f) N(\tau,t) - R_{12}^{(n)}(f) H_{2}^{s}(f) S(\tau,f) \\
\\
= \ &  \underbrace{S(\tau,f) \left( H_{1}^{s}(f) - R_{12}^{(n)}(f) H_{2}^{s}(f) \right)}_{= \hat{S}(\tau,f)}
\end{split}
$$
which shows that the noise's signal $N(\tau,f)$ gets cancelled off and a filtered speech signal remains. 

In this sense, we can denoise the target microphone's recording by subtracting an estimate of the noise component obtained from the reference microphone and known ReTF, denoted by
$$
\begin{equation}
\hat{S}(\tau,f) \approx P_{1}(\tau,f) - R_{12}^{(n)}(f) P_{2}(\tau,f)
\end{equation}
$$
where $\hat{S}(\tau,f)$ can be considered an ehancement (with slight distortion) of whatever signal other than the noise is present in the target microphone, in this case the speech signal.

### Estimating the ReTF for Denoising
The signal denoising of (3) requires the ReTF of the noise source ($n$) to be known. In practice, we can estimate the noise's ReTF by the assumption that on average over time the recording contains mostly active noise and inactive speech. Therefore, the spectral covariance approach in [Obtaining ReTF in Practice](#obtaining-retf-in-practice) can be used such that 
$$
R_{1,2}^{(n)}(f) 
\approx \frac{ \mathcal{P}_{1,2}(f) }{ \mathcal{P}_{2,2}(f) }
\approx \frac{ \frac{1}{T} \sum_{t=1}^{T} \Rho_{1,2}(\tau-t,f) }{ \frac{1}{T} \sum_{t=1}^{T} \Rho_{2,2}(\tau-t,f) } 
% \approx \frac{ H_{1}^{n}(f) \tilde{N} H_{2}^{n}(f)^{*} }{ H_{2}^{n}(f) \tilde{N} H_{2}^{n}(f)^{*} } \approx \frac{ H_{1}^{n}(f) }{ H_{2}^{n}(f) }
$$
because
$$
\frac{ \frac{1}{T} \sum_{t=1}^{T} \Rho_{1,2}(\tau-t,f) }{ \frac{1}{T} \sum_{t=1}^{T} \Rho_{2,2}(\tau-t,f) } \approx \frac{ H_{1}^{n}(f) \tilde{N} H_{2}^{n}(f)^{*} }{ H_{2}^{n}(f) \tilde{N} H_{2}^{n}(f)^{*} } \approx \frac{ H_{1}^{n}(f) }{ H_{2}^{n}(f) }
$$
where $\tilde{N}$ denotes the PSD of the noise.

### Step-by-step ReTF Denoising

1. Load / obtain the recordings $p_{1}(t)$ and $p_{2}(t)$.

2. STFT the recordings to get $P_{1}(\tau,f)$ and $P_{2}(\tau,f)$. Use a window length greater than the RT60 of the room to insure the window length is > $h_{\{1,2\}}^{n}(t)$. Because we are filtering in the STFT domain in (3), we also need to use a zero-padded NFFT size greater than 2*window length to avoid circular convolution. Use a hop size of window length / 2 or better.
3. Estimate the CPSDs of $\Rho_{1,2}(\tau,f)$ and $\Rho_{2,2}(\tau,f)$.
4. Get the average CPSD (mean) over the whole recording $\mathcal{P}_{1,2}(f)$ and $\mathcal{P}_{2,2}(f)$. Because we assume that the noise is active and speech is inactive on averge, these averaged CPSDs are just of the noise source $\mathcal{P}_{1,2}(f) \approx \mathcal{P}_{1,2}^{(n)}(f)$ and $\mathcal{P}_{2,2}(f) \approx \mathcal{P}_{2,2}^{(n)}(f)$.
5. Estimate the noise's ReTF: 
$$
R_{1,2}^{(n)}(f) 
\approx \frac{ \mathcal{P}_{1,2}^{(n)}(f) }{ \mathcal{P}_{2,2}^{(n)}(f) }
$$
6. Filter the reference microphone's recording with the ReTF to estimate the noise component of the target microphone:
$$
P_{1}^{(n)}(\tau,f) \approx R_{12}^{(n)}(f) \times P_{2}(\tau,f)
$$
7. Subtract the estimated noise component from the target microphone's recording to obtain a denoised the speech estimate $\hat{S}(\tau,f)$:
$$
\hat{S}(\tau,f) \approx P_{1}(\tau,f) - P_{1}^{(n)}(\tau,f)
$$
8. ISTFT overlap-add $\hat{S}(\tau,f)$ to return to time-domain signal of denoised speech $\hat{s}(t)$. Note that the ISTFT can also be applied before step-7 to do the noise subtraction in the time-domain:
$$
\hat{s}(t) \approx p_{1}(t) - p_{1}^{(n)}(t)
$$

















## License

This project is licensed under the [BSD 3-Clause License](LICENSE).
