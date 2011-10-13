-- | BAM - extract information from BAM files

import Classify
import Bio.SamTools.Bam.Unsafe
import System.Environment (getArgs)

main = do
  [f] <- getArgs
  putStrLn . display . classify . take 10000 =<< readBams f