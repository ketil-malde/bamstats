-- | BAM - extract information from BAM files

{-# Language DeriveDataTypeable #-}

module Main where
import Classify
import Bio.SamTools.Bam.Unsafe
import System.Console.CmdArgs

data Options = Stats { numrd :: Maybe Int, inputs :: [FilePath] }
             | Hist { numrd :: Maybe Int, bins, maxdist :: Int, inputs :: [FilePath] }
             deriving (Data,Typeable,Read,Show)

clfy :: Options
clfy = Stats { numrd = Nothing &= help "max number of reads to include"
                , inputs = []     &= args &= typ "BAM file(s)"
                } &= help "Calculate statistics on insert sizes"
hst  = Hist { numrd = Nothing &= help "max number of reads to include"
            , bins  = 20      &= help "number of bins" &= typ "INT"
            , maxdist = 1000  &= help "max insert size to count" &= typ "INT"
            , inputs = []     &= args &= typ "BAM file(s)"
            } &= help "Collect a histogram of insert sizes"

main = do
  o <- cmdArgs $ modes [clfy,hst] &= help "Extract information from BAM files" 
       &= program "bam" &= summary "bam v0.0, ©2011 Ketil Malde"
  let geninp f = (case numrd o of Just x -> take x; Nothing -> id) `fmap` readBams f
      genout = case o of Stats {} -> putStrLn . display . classify
                         Hist {}     -> const (print o)
  mapM_ (\f -> genout =<< geninp f) $ inputs o