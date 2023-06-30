# QKD
Python scripts for simulating quantum key distribution protocols (primarily BB84) and attacks used against it (intercept and resend, photon-number-splitting, Slutsky-Brandt attack). These scripts were written in 2019-2021 by Leo Crowder as an undergraduate at Northern Arizona University as a research project under the direction of Ines Montanto. 

The BB84 Protocol, originially proposed by Charles Bennett and Gilles Brassard in 1984, was the first quantum cryptography protocol https://en.wikipedia.org/wiki/BB84. BB84 outlines a procedure for two parties, Alice and Bob, to use quantum mechanics to establish a private key for one-time pad encryption. The fundamental idea is that any attempts made to intercept their communications will necessarily require altering the quantum states of the signal, which results in a highly noticeable disturbance (discernable from the expected errors due to the randomness of quantum measurements).

In practice, QKD is typically implemented using photons, and have vulnerabilities that are not accounted for by the theory. There are various errors in equipment such as channel loss, channel "drift" (where the quantum state is altered), and measurement noise/errors. This opens up new opportunities for eavesdroppers (Eve) to exploit these imperfections. Much of the code here is an exploration of attack methods (intercept and resend, photon-number splitting, Slutsky-Brandt) and general eavesdropping strategies to maximize information gain while going undetected. While not particularly advanced due to my experience level at the time, there is some statistical analysis of these attacks (hypothesis testing, power analysis). 

I started this project with minimal Python experience and learned a significant amount along the way, so the code will have varying degrees of quality. The most recent and most competent work within this project is in the `netsquid` directory. In the summer of 2021 I used NetSquid (https://netsquid.org), a quantum network design software, to implement a decoy-state protocol with realistic quantum hardware components. See the `writeups` a summary of this work, as well as a poster for my 2021 Hooper Undergraduate Research Award presentation.
