void chain_justinGen_pythia(Int_t nev = 10000, char* filename = "outloop_io_test_10k.root", char* folderName){

    // Runloader
    TStopwatch timer;
    timer.Start();

    gRandom->SetSeed(0);//put 0 to use system time
    gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT");

    gSystem->Load("libhijing.so");
    gSystem->Load("libTHijing.so");     // AliGenHijing is defined here
    gSystem->Load("liblhapdf.so");      // Parton density functions
    gSystem->Load("libEGPythia6.so");   // TGenerator interface
    gSystem->Load("libpythia6.so");     // Pythia
    gSystem->Load("libAliPythia6.so");  // ALICE specific implementations
    gSystem->Load("libEVGEN.so");      // EVGEN library to allow for "SetTriggerParticle" function


    AliRunLoader* rl = AliRunLoader::Open(filename, folderName, "recreate");


    string kineFileName = "Kinematics_";
    kineFileName.append(filename);
    rl->SetKineFileName(kineFileName);
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(nev);
    rl->LoadKinematics("RECREATE");
    rl->MakeTree("E");

    gAlice->SetRunLoader(rl);

    //Create stack
    rl->MakeStack();
    AliStack* stack = rl->Stack();

    //Header
    AliHeader* header = rl->GetHeader();

    AliGenPythia * gener = new AliGenPythia(-1);   something like event scaling
    gener->SetEnergyCMS(7000.);
//    gener->SetNuclei(208,208);
//    gener->SetPtHard(10,10000);
    gener->SetOrigin(0.,0.,0.); 
    gener->SetSigma(0.,0.,5.3);
    gener->SetVertexSmear(kPerEvent);
//    gener->SetTriggerParticle(0, 1., 5., 1000.);
//    gener->SetChildPtRange(5, 1000);

    //Initialize generator
    gener->Init();
    gener->SetStack(stack);

    //Event Loop

    Int_t iev;
    for (iev = 0; iev < nev; iev++) {
        if (iev % 50 == 0)  Printf(Form(" Event number %d", iev));
        //Initialize event
        header->Reset(0,iev);
        rl->SetEventNumber(iev);
        stack->Reset();
        rl->MakeTree("K");
        rl->MakeStack();
        stack->ConnectTree(rl->TreeK());

        //Generate event
        gener->Generate();

        //Analysis
        Int_t npart = stack->GetNprimary();

        //Finish event
        header->SetNprimary(stack->GetNprimary());
        header->SetNtrack(stack->GetNtrack());

        //I/O
        stack->FinishEvent();
        header->SetStack(stack);
        rl->TreeE()->Fill();
        rl->WriteKinematics("OVERWRITE");

    }   //end event loop

    //Termination

    //Generator
    gener->FinishRun();

    //Write file
    rl->WriteHeader("OVERWRITE");
    gener->Write();
    rl->Write();

    rl->RemoveEventFolder();
    rl->UnloadAll();

    timer.Stop();
    timer.Print();
}
