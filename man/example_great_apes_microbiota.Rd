\name{example_great_apes_microbiota}
\alias{example_great_apes_microbiota}

\title{
Example of empirical application of HOME
}

\description{
This function provides an example of empirical application using HOME for 3 different OTUs from the great apes microbiota. It only downloads the data (alignments and host tree).
}

\usage{example_great_apes_microbiota(name,path)}

\arguments{
\item{name}{
the name of the run
}

\item{path}{
optional - the path to the working directory
}
}

\references{
Perez-Lamarque, B. and Morlon, H. (in. prep.) Characterizing symbiont inheritance during host-microbiota evolution : application to the great apes microbiota
}
\author{Benoît Perez-Lamarque}
\seealso{
\code{\link{Home_model}}
}
\examples{
# This example may take a little bit of time. Be patient!

example_great_apes_microbiota(name="great_apes")

# Run HOME (parameters estimations and model selections)

#HOME_model(name="great_apes", name_index=c("OTU0001","OTU0002","OTU0003"), nb_tree=500, lambda=seq(1,25), nb_random=10)

# NB: this should be run on a larger number of tree!

}